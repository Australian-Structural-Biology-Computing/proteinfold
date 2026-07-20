ARG CUDA=12.6.0
FROM nvidia/cuda:${CUDA}-base-ubuntu22.04

# FROM directive resets ARGS, so we specify again (the value is retained if
# previously set).
ARG CUDA

LABEL org.opencontainers.image.title="nf-core/proteinfold_chai1" \
    org.opencontainers.image.description="Docker image containing Chai-1 for nf-core/proteinfold" \
    org.opencontainers.image.version="dev" \
    org.opencontainers.image.authors="Joshua Caley <j.caley@unsw.edu.au>" \
    org.opencontainers.image.source="https://github.com/nf-core/proteinfold" \
    org.opencontainers.image.licenses="MIT"

ENV VIRTUAL_ENV=/opt/venv
ENV PATH="${VIRTUAL_ENV}/bin:${PATH}"

RUN apt-get update --quiet \
    && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends --yes --quiet \
        kalign \
        python3 \
        python3-pip \
        python3-venv \
    && python3 -m venv "${VIRTUAL_ENV}" \
    && pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir chai_lab==0.6.1 \
    && rm -rf /var/lib/apt/lists/*
