#!/usr/bin/env python3
"""Cluster predicted structures of one sample and project them onto a 2D plot.

All ranked structures predicted for a single sample (across every prediction
backend, model rank and seed) are compared pairwise, grouped into
conformational clusters and embedded in 2D with classical multidimensional
scaling (MDS).

Structural distances come from a Foldseek-Multimer all-vs-all complex
alignment when available: Foldseek optimises the chain pairing between two
complexes, so the complex TM-scores are invariant to chain naming, ordering
and composition differences between prediction backends. Without Foldseek
results an approximate TM-score is computed from a Kabsch superposition of
the CA atoms instead (chain correspondence is then taken from the order of
appearance and residue numbering).

Outputs:
  * `<name>_conformations.tsv`   per-structure table (model, rank, seed,
                                 cluster, confidence metrics, 2D coordinates)
  * `<name>_conformations.html`  interactive 2D scatter plot with options to
                                 color and filter by model and confidence score

Clustering joins two structures when their complex TM-score similarity is at
or above `--cluster-cutoff`.
"""

import argparse
import csv
import json
import os
import re
import sys

import numpy as np

from Bio.PDB import MMCIFParser, PDBParser

SUFFIX_RE = re.compile(r"rank_(\d+).*?seed_(\d+)(?:_sample_(\d+))?")
METRICS = ["ptm", "iptm", "ipsae", "ranking_score", "plddt"]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structures", nargs="+", required=True,
                        help="Predicted structure files (PDB or mmCIF) for one sample.")
    parser.add_argument("--scores", nargs="*", default=[],
                        help="Per-structure ColabFold confidence JSON files.")
    parser.add_argument("--manifest", required=True,
                        help="TSV mapping each structure file to its prediction model: name, model, mode_short.")
    parser.add_argument("--multimer-report",
                        help="Foldseek-Multimer all-vs-all complex report for the sample (query, target, matched chains, TM-scores, ...).")
    parser.add_argument("--name", required=True, help="Sample identifier (meta.id).")
    parser.add_argument("--cluster-cutoff", type=float, default=0.5,
                        help="Complex TM-score similarity at or above which two structures are joined into the same conformational cluster.")
    parser.add_argument("--output-dir", default=".")
    return parser.parse_args()


def structure_suffix(name):
    """Suffix shared by a structure file and its scores JSON (`rank_..._sample_N`).

    ColabFold names raw outputs after the FASTA query name, which can differ
    from the sample id, so the suffix is taken from the last `_rank_` token
    rather than from a name prefix.
    """
    stem = os.path.splitext(name)[0]
    head, sep, tail = stem.rpartition("_rank_")
    if not sep or not re.match(r"\d", tail):
        return None
    return f"rank_{tail}"


def load_manifest(path, structure_names):
    """Map structure file name -> (model, mode_short) from the manifest TSV."""
    labels = {}
    with open(path) as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) >= 2:
                labels[row[0]] = (row[1], row[2] if len(row) > 2 else row[1])
    return {name: labels.get(name, ("unknown", "unknown")) for name in structure_names}


def ca_coordinates(path):
    """CA atom coordinates keyed by (chain order of appearance, residue number)."""
    ext = os.path.splitext(path)[1].lower()
    if ext == ".pdb":
        parser = PDBParser(QUIET=True)
    elif ext == ".cif":
        parser = MMCIFParser(QUIET=True)
    else:
        raise ValueError(f"Unsupported structure format: {path}")
    structure = parser.get_structure(os.path.basename(path), path)
    model = next(iter(structure))
    coords, b_factors = {}, {}
    for chain_idx, chain in enumerate(model):
        for residue in chain:
            if residue.id[0] != " ":  # skip HETATM/water records
                continue
            ca = residue["CA"] if "CA" in residue else None
            if ca is None:
                continue
            key = (chain_idx, residue.id[1])
            coords[key] = ca.coord.astype(float)
            b_factors[key] = float(ca.bfactor)
    return coords, b_factors


def read_metrics(scores_by_suffix, suffix, structure_bfactors):
    """Confidence metrics for one structure from its scores JSON (may be absent)."""
    metrics = {metric: None for metric in METRICS}
    plddt_values = [value for value in structure_bfactors.values() if value > 0]
    if plddt_values:
        metrics["plddt"] = float(np.mean(plddt_values))
    scores = scores_by_suffix.get(suffix)
    if not scores:
        return metrics
    for metric in ("ptm", "iptm", "ranking_score"):
        value = scores.get(metric)
        if isinstance(value, (int, float)):
            metrics[metric] = float(value)
    ipsae = scores.get("ipsae")
    if isinstance(ipsae, dict) and ipsae:
        metrics["ipsae"] = float(max(ipsae.values()))
    elif isinstance(ipsae, (int, float)):
        metrics["ipsae"] = float(ipsae)
    if metrics["plddt"] is None:
        plddt = scores.get("plddt")
        if isinstance(plddt, (list, tuple)) and plddt:
            metrics["plddt"] = float(np.mean(plddt))
        elif isinstance(plddt, (int, float)):
            metrics["plddt"] = float(plddt)
    return metrics


def read_multimer_report(path, names):
    """Parse a Foldseek-Multimer complex report into a similarity map.

    Report columns: query complex, target complex, matched query chains,
    matched target chains, TM-score normalised by the query length, TM-score
    normalised by the target length, rotation matrix, translation vector and
    the assignment id. A pair may carry several alternative chain assignments;
    the best (highest) score is kept and the similarity is the smaller of the
    two length-normalised TM-scores, matching Foldseek's MMscore convention.
    Query and target identifiers are structure file names without extension.
    """
    stems = {os.path.splitext(name)[0]: name for name in names}
    similarity, assignments = {}, {}
    with open(path) as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 6 or row[0] not in stems or row[1] not in stems:
                continue
            try:
                score = min(float(row[4]), float(row[5]))
            except ValueError:
                continue
            pair = (row[0], row[1])
            if pair not in similarity or score > similarity[pair]:
                similarity[pair] = score
                assignments[pair] = (row[2], row[3])
    return similarity, assignments


def kabsch_rmsd(coords_a, coords_b):
    """RMSD and approximate TM-score after optimal superposition.

    The TM-score is estimated from the superposed CA distances with the
    standard TM-score denominator d0; the rotation is RMSD-optimal rather
    than TM-optimal, so the value is approximate.
    """
    center_a = coords_a.mean(axis=0)
    center_b = coords_b.mean(axis=0)
    pa, pb = coords_a - center_a, coords_b - center_b
    u, _, vt = np.linalg.svd(pa.T @ pb)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rotation = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    diff = (rotation @ pa.T).T - pb
    rmsd = float(np.sqrt((diff ** 2).sum(axis=1).mean()))
    length = len(pa)
    d0 = 1.24 * (max(length, 20) - 15) ** (1 / 3) - 1.8
    tm_score = float(np.sum(1.0 / (1.0 + (np.linalg.norm(diff, axis=1) / d0) ** 2)) / length)
    return rmsd, tm_score


def pairwise_distances(uploads, similarity=None):
    """Symmetric distance matrix (1 - complex TM-score) over all structures.

    With a Foldseek-Multimer similarity map, pairs without an alignment are
    treated as maximally distant. Without one, an approximate TM-score is
    computed from a Kabsch superposition of the common CA atoms (chain
    correspondence from order of appearance and residue numbering).
    """
    count = len(uploads)
    matrix = np.full((count, count), 1.0)
    np.fill_diagonal(matrix, 0.0)
    stems = [os.path.splitext(uploads[i][0])[0] for i in range(count)]
    for i in range(count):
        for j in range(i + 1, count):
            if similarity is not None:
                pair_score = similarity.get((stems[i], stems[j]), similarity.get((stems[j], stems[i])))
                if pair_score is None:
                    continue  # no alignment found: keep the maximal distance
            else:
                shared = sorted(set(uploads[i][1]) & set(uploads[j][1]))
                if len(shared) < 3:
                    continue  # no meaningful comparison: keep the maximal distance
                a = np.array([uploads[i][1][key] for key in shared])
                b = np.array([uploads[j][1][key] for key in shared])
                _, pair_score = kabsch_rmsd(a, b)
            matrix[i, j] = matrix[j, i] = 1.0 - min(max(pair_score, 0.0), 1.0)
    return matrix


def average_linkage_clusters(matrix, max_distance):
    """Greedy average-linkage agglomerative clustering cut at `max_distance`."""
    n = matrix.shape[0]
    clusters = [[i] for i in range(n)]
    while len(clusters) > 1:
        best, best_pair = None, None
        for a in range(len(clusters)):
            for b in range(a + 1, len(clusters)):
                dist = matrix[np.ix_(clusters[a], clusters[b])].mean()
                if best is None or dist < best:
                    best, best_pair = dist, (a, b)
        if best is None or best > max_distance:
            break
        a, b = best_pair
        clusters[a] = clusters[a] + clusters[b]
        del clusters[b]
    # Name clusters by descending size: cluster_1 is the largest.
    order = sorted(range(len(clusters)), key=lambda idx: (-len(clusters[idx]), clusters[idx][0]))
    labels = {}
    for rank, idx in enumerate(order, start=1):
        for member in clusters[idx]:
            labels[member] = f"cluster_{rank}"
    return labels


def classical_mds(matrix):
    """Classical MDS embedding of a distance matrix into 2D."""
    n = matrix.shape[0]
    if n == 1:
        return np.zeros((1, 2))
    squared = matrix ** 2
    centering = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * centering @ squared @ centering
    b = (b + b.T) / 2  # symmetrise numerical noise
    eigenvalues, eigenvectors = np.linalg.eigh(b)
    order = np.argsort(eigenvalues)[::-1][:2]
    values = np.clip(eigenvalues[order], 0.0, None)
    coords = eigenvectors[:, order] * np.sqrt(values)
    if coords.shape[1] < 2:
        coords = np.hstack([coords, np.zeros((n, 2 - coords.shape[1]))])
    return coords


def build_html(sample, cutoff, points, path, used_foldseek=True):
    """Self-contained interactive 2D conformational plot (plotly via CDN)."""
    context = {
        "sample": sample,
        "cutoff": cutoff,
        "metrics": METRICS,
        "points": points,
    }
    payload = json.dumps(context).replace("</", "<\\/")
    html = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Conformational analysis - __SAMPLE__</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
  body { font-family: sans-serif; margin: 0; }
  header { padding: 12px 20px; border-bottom: 1px solid #ddd; }
  header h1 { font-size: 18px; margin: 0 0 4px 0; }
  header p { margin: 0; color: #555; font-size: 13px; }
  #controls { display: flex; gap: 24px; flex-wrap: wrap; align-items: center;
              padding: 10px 20px; border-bottom: 1px solid #ddd; }
  #controls label { font-size: 13px; color: #333; }
  #controls select { margin-left: 6px; }
  #filter-value-wrap { display: none; align-items: center; gap: 8px; }
  #plot { width: 100%; height: 82vh; }
</style>
</head>
<body>
<header>
  <h1>Conformational analysis: __SAMPLE__</h1>
  <p>__N__ predicted structures clustered by __METHOD__ (average linkage, similarity &ge; __CUTOFF__) and embedded in 2D with classical MDS.</p>
</header>
<div id="controls">
  <label>Color by
    <select id="color-by"></select>
  </label>
  <label>Filter by
    <select id="filter-by">
      <option value="">none</option>
    </select>
  </label>
  <span id="filter-value-wrap">
    <label>min <span id="filter-metric"></span></label>
    <input type="range" id="filter-slider" min="0" max="1" step="0.01" value="0">
    <span id="filter-value"></span>
  </span>
</div>
<div id="plot"></div>
<script type="application/json" id="conformations-data">__PAYLOAD__</script>
<script>
const data = JSON.parse(document.getElementById('conformations-data').textContent);
const CATEGORICAL = ['cluster', 'model', 'rank', 'seed'];
const PALETTES = {};
const QUALITATIVE = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
                     '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
                     '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5'];
function palette(field) {
  if (!PALETTES[field]) {
    const values = [...new Set(data.points.map(p => String(p[field])))].sort();
    PALETTES[field] = Object.fromEntries(values.map((v, i) => [v, QUALITATIVE[i % QUALITATIVE.length]]));
  }
  return PALETTES[field];
}
const colorSelect = document.getElementById('color-by');
const filterSelect = document.getElementById('filter-by');
const filterWrap = document.getElementById('filter-value-wrap');
const filterSlider = document.getElementById('filter-slider');
const filterValue = document.getElementById('filter-value');
for (const field of CATEGORICAL.concat(data.metrics)) {
  colorSelect.add(new Option(field, field));
  if (CATEGORICAL.includes(field)) continue;
  filterSelect.add(new Option(field, field));
}
function fmt(value) { return value == null ? 'n/a' : value.toFixed(3); }
function visiblePoints() {
  const metric = filterSelect.value;
  if (!metric) return data.points;
  const threshold = parseFloat(filterSlider.value);
  return data.points.filter(p => p.metrics[metric] != null && p.metrics[metric] >= threshold);
}
function buildFigure() {
  const field = colorSelect.value;
  const points = visiblePoints();
  const hover = points.map(p =>
    `<b>${p.name}</b><br>model: ${p.model}<br>rank: ${p.rank}, seed: ${p.seed}` +
    data.metrics.map(m => `<br>${m}: ${fmt(p.metrics[m])}`).join(''));
  let marker, traces;
  if (CATEGORICAL.includes(field)) {
    const colors = palette(field);
    const groups = {};
    points.forEach((p, i) => {
      const key = String(p[field]);
      (groups[key] = groups[key] || {x: [], y: [], text: []});
      groups[key].x.push(p.x); groups[key].y.push(p.y); groups[key].text.push(hover[i]);
    });
    traces = Object.keys(groups).sort().map(key => ({
      type: 'scattergl', mode: 'markers', name: key,
      x: groups[key].x, y: groups[key].y, text: groups[key].text,
      hoverinfo: 'text',
      marker: {size: 9, color: colors[key], line: {width: 1, color: '#333'}},
    }));
    marker = null;
  } else {
    const values = points.map(p => p.metrics[field] == null ? null : p.metrics[field]);
    const finite = values.filter(v => v != null);
    const zmin = finite.length ? Math.min(...finite) : 0;
    const zmax = finite.length ? Math.max(...finite) : 1;
    traces = [{
      type: 'scattergl', mode: 'markers',
      x: points.map(p => p.x), y: points.map(p => p.y), text: hover,
      hoverinfo: 'text',
      marker: {size: 9, color: values, colorscale: 'Viridis', cmin: zmin, cmax: zmax,
               line: {width: 1, color: '#333'},
               colorbar: {title: field}},
    }];
    marker = {showlegend: false};
  }
  return {
    data: traces,
    layout: {
      margin: {t: 10, b: 40, l: 40, r: 10},
      xaxis: {title: 'MDS 1'}, yaxis: {title: 'MDS 2', scaleanchor: 'x'},
      hovermode: 'closest',
      showlegend: marker === null ? true : false,
      legend: {itemsizing: 'constant'},
    },
    config: {responsive: true, displayModeBar: false},
  };
}
function render() { Plotly.react('plot', buildFigure().data, buildFigure().layout, buildFigure().config); }
colorSelect.addEventListener('change', render);
filterSelect.addEventListener('change', () => {
  const metric = filterSelect.value;
  filterWrap.style.display = metric ? 'inline-flex' : 'none';
  if (!metric) { render(); return; }
  const values = data.points.map(p => p.metrics[metric]).filter(v => v != null);
  const lo = values.length ? Math.min(...values) : 0;
  const hi = values.length ? Math.max(...values) : 1;
  filterSlider.min = lo; filterSlider.max = hi;
  filterSlider.step = (hi - lo) / 100 || 0.01;
  filterSlider.value = lo;
  document.getElementById('filter-metric').textContent = metric;
  updateFilter();
});
function updateFilter() {
  filterValue.textContent = parseFloat(filterSlider.value).toFixed(3);
  render();
}
filterSlider.addEventListener('input', updateFilter);
render();
</script>
</body>
</html>"""
    html = (html
            .replace("__SAMPLE__", sample)
            .replace("__N__", str(len(points)))
            .replace("__METHOD__", "Foldseek-Multimer complex TM-score" if used_foldseek
                     else "approximate TM-score from CA superposition")
            .replace("__CUTOFF__", f"{cutoff:g}")
            .replace("__PAYLOAD__", payload))
    with open(path, "w") as handle:
        handle.write(html)


def main():
    args = parse_args()
    structure_names = [os.path.basename(path) for path in args.structures]
    labels = load_manifest(args.manifest, structure_names)

    scores_by_suffix = {}
    for path in args.scores:
        name = os.path.basename(path)
        match = re.search(r"_scores_(rank_\d+.*)\.json$", name)
        if match:
            try:
                with open(path) as handle:
                    scores_by_suffix[match.group(1)] = json.load(handle)
            except (OSError, json.JSONDecodeError) as error:
                print(f"WARNING: could not read scores from {name}: {error}", file=sys.stderr)

    rows, uploads = [], []
    for path in args.structures:
        name = os.path.basename(path)
        suffix = structure_suffix(name)
        rank = seed = sample = None
        if suffix:
            match = SUFFIX_RE.search(suffix)
            if match:
                rank, seed, sample = int(match.group(1)), int(match.group(2)), int(match.group(3) or 0)
        try:
            coords, b_factors = ca_coordinates(path)
        except Exception as error:
            print(f"WARNING: could not parse structure {name}: {error}", file=sys.stderr)
            continue
        if not coords:
            print(f"WARNING: no CA atoms found in {name}", file=sys.stderr)
            continue
        model, mode_short = labels[name]
        rows.append({
            "structure": name,
            "model": model,
            "mode_short": mode_short,
            "rank": rank,
            "seed": seed,
            "sample": sample,
            **read_metrics(scores_by_suffix, suffix, b_factors),
        })
        uploads.append((name, coords))
    if not uploads:
        sys.exit("CLUSTER_CONFORMATIONS found no parseable predicted structures")

    similarity = {}
    if args.multimer_report:
        try:
            similarity, _ = read_multimer_report(args.multimer_report, [u[0] for u in uploads])
        except OSError as error:
            print(f"WARNING: could not read the Foldseek-Multimer report: {error}", file=sys.stderr)
    if not similarity:
        print("WARNING: no Foldseek-Multimer similarities available; "
              "falling back to approximate TM-scores from a Kabsch superposition", file=sys.stderr)
    matrix = pairwise_distances(uploads, similarity or None)
    cluster_labels = average_linkage_clusters(matrix, 1.0 - args.cluster_cutoff)
    coords_2d = classical_mds(matrix)

    for index, row in enumerate(rows):
        row["cluster"] = cluster_labels[index]
        row["mds1"] = round(float(coords_2d[index, 0]), 3)
        row["mds2"] = round(float(coords_2d[index, 1]), 3)

    os.makedirs(args.output_dir, exist_ok=True)
    tsv_path = os.path.join(args.output_dir, f"{args.name}_conformations.tsv")
    html_path = os.path.join(args.output_dir, f"{args.name}_conformations.html")

    columns = ["structure", "model", "mode_short", "rank", "seed", "sample", "cluster",
               *METRICS, "mds1", "mds2"]
    with open(tsv_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t",
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: ("" if row[key] is None else row[key]) for key in columns})

    points = [{
        "name": row["structure"],
        "model": row["model"],
        "rank": row["rank"],
        "seed": row["seed"],
        "cluster": row["cluster"],
        "x": row["mds1"],
        "y": row["mds2"],
        "metrics": {metric: row[metric] for metric in METRICS},
    } for row in rows]
    build_html(args.name, args.cluster_cutoff, points, html_path, bool(similarity))


if __name__ == "__main__":
    main()
