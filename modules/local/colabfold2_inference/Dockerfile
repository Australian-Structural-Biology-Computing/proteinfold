ARG CUDA_VERSION=12.9.1

# MMseqs2, so the same image can also drive colabfold_search (the MSA search
# processes share this runtime). Built from the pipeline's MMseqs2 fork, whose
# pairing behaviour the search relies on, matching the mmseqs_colabfoldsearch
# image.
FROM --platform=linux/amd64 ghcr.io/steineggerlab/build-containers:main-x86_64 AS mmseqs-builder

ARG MMSEQS2_REPOSITORY=https://github.com/jscgh/MMseqs2.git
ARG MMSEQS2_REF=48f2a22d5f68f595afdfb5e601ca11721cf17449

WORKDIR /opt/mmseqs2

RUN git init . && \
    git remote add origin "${MMSEQS2_REPOSITORY}" && \
    git fetch --depth 1 origin "${MMSEQS2_REF}" && \
    git checkout --detach FETCH_HEAD && \
    cmake -S . -B build -GNinja \
        -DCMAKE_BUILD_TYPE=Release \
        -DHAVE_TESTS=0 \
        -DHAVE_SSE4_1=1 \
        -DFORCE_STATIC_DEPS=1 && \
    cmake --build build -j"$(nproc)"

FROM nvidia/cuda:${CUDA_VERSION}-runtime-ubuntu24.04

ARG COLABFOLD_COMMIT=38683c21700408b6a29cf395e8ef07cf109fcb91

LABEL org.opencontainers.image.title="nf-core/proteinfold_colabfold2" \
      org.opencontainers.image.description="ColabFold AF3 preview runtime for ProteinFold" \
      org.opencontainers.image.version="af3-preview-38683c2" \
      org.opencontainers.image.source="https://github.com/sokrypton/ColabFold/tree/38683c21700408b6a29cf395e8ef07cf109fcb91" \
      org.opencontainers.image.revision="38683c21700408b6a29cf395e8ef07cf109fcb91" \
      org.opencontainers.image.licenses="MIT; model weights have separate terms"

COPY requirements.lock /opt/colabfold-requirements.lock

# Keep the runtime pins fixed so rebuilding this source cannot silently pick a different JAX runtime or transitive dependency.
RUN apt-get update && \
    apt-get install -y --no-install-recommends git python3 python3-pip python3-venv && \
    python3 -m venv /opt/colabfold && \
    /opt/colabfold/bin/pip install --no-cache-dir "pip==26.2.1" && \
    /opt/colabfold/bin/pip install --no-cache-dir -r /opt/colabfold-requirements.lock && \
    apt-get purge -y git && \
    apt-get autoremove -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /root/.cache

COPY --from=mmseqs-builder /opt/mmseqs2/build/src/mmseqs /opt/colabfold/bin/mmseqs

RUN /opt/colabfold/bin/mmseqs version

# Prebuild alphafold3's chemical component tables, which otherwise get written into site-packages
# on first use and fail under read-only runtimes; components.cif is staged via params/af3/libcifpp/.
RUN mkdir /tmp/ccd_build && \
    /opt/colabfold/bin/python -c "from colabfold.alphafold3 import require; require('/tmp/ccd_build'); require('/tmp/ccd_build')" && \
    rm -rf /tmp/ccd_build

ENV PATH="/opt/colabfold/bin:${PATH}" \
    MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/mplconfig \
    XDG_CACHE_HOME=/tmp/xdg_cache \
    JAX_COMPILATION_CACHE_DIR=/tmp/jax_cache
