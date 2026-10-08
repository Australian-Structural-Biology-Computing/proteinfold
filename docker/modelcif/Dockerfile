FROM python:3.12-slim

RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    "biopython==1.84" \
    "pyyaml==6.0.2" \
    "modelcif==1.7" \
    "msgpack>=1.0"
