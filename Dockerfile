FROM mambaorg/micromamba:1.5.10

COPY environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
    && micromamba clean --all --yes

ENV PATH=/opt/conda/bin:$PATH
ENV MAMBA_DOCKERFILE_ACTIVATE=1

# Nextflow collects task metrics with `ps`.
RUN micromamba install -y -n base -c conda-forge procps-ng \
    && micromamba clean --all --yes

# Nextflow needs a working bash; the micromamba base image provides it.
SHELL ["/bin/bash", "-lc"]
WORKDIR /work

# Nextflow mounts the work directory from the host. Running as root and using
# docker.fixOwnership = true avoids host/container UID mismatches.
USER root
