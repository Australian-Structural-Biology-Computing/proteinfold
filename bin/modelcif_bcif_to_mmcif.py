#!/usr/bin/env python3
"""Convert BinaryCIF to text ModelCIF for dictionary validation."""

import sys

import modelcif.dumper
import modelcif.reader

with open(sys.argv[1], "rb") as source:
    systems = modelcif.reader.read(source, format="BCIF")
with open(sys.argv[2], "w") as destination:
    modelcif.dumper.write(destination, systems, format="mmCIF")
