FROM mambaorg/micromamba:1.5.10

COPY environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
    && micromamba clean --all --yes

ENV PATH=/opt/conda/bin:$PATH
ENV MAMBA_DOCKERFILE_ACTIVATE=1

# Nextflow needs a working bash; the micromamba base image provides it.
SHELL ["/bin/bash", "-lc"]
WORKDIR /work
