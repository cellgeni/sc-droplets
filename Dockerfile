# Base image: Python 3.12 + uv preinstalled (Debian slim), pinned by digest
# (of the multi-arch index) so a rebuilt release gets the same base. To move
# it, take the index digest of the new tag:
#   docker buildx imagetools inspect ghcr.io/astral-sh/uv:<tag>
FROM ghcr.io/astral-sh/uv:0.9.30-python3.12-bookworm-slim@sha256:e5b65587bce7de595f299855d7385fe7fca39b8a74baa261ba1b7147afa78e58

# PYTHONNOUSERSITE: Apptainer bind-mounts the host $HOME by default, so a user's
#   ~/.local/lib/python3.12/site-packages would otherwise shadow this venv.
# XDG_CACHE_HOME, MPLCONFIGDIR: Nextflow is often run with `-u $(id -u):$(id -g)`,
#   which leaves the container with no writable $HOME; matplotlib needs a
#   writable config dir.
# MPLBACKEND: no display in a container.
ENV UV_NO_DEV=1 \
    UV_COMPILE_BYTECODE=1 \
    PYTHONNOUSERSITE=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    XDG_CACHE_HOME=/tmp/.cache \
    MPLCONFIGDIR=/tmp/.cache/matplotlib \
    MPLBACKEND=Agg

# procps supplies `ps`, which Nextflow needs to collect per-task metrics.
RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates procps \
    && rm -rf /var/lib/apt/lists/*

# Fail the build, rather than every Nextflow task, if the base image ever drops
# one of the tools Nextflow requires in a task container.
RUN set -eu; for t in bash ps awk date grep sed tail tee; do \
      command -v "$t" >/dev/null || { echo "missing required tool: $t" >&2; exit 1; }; \
    done

WORKDIR /cli

COPY . .

# --locked asserts that uv.lock is in sync with pyproject.toml. The cache is
# cleared in the same layer and replaced by an empty world-writable folder,
# so a task running under any UID can write to the path XDG_CACHE_HOME names.
RUN uv sync --locked \
    && rm -rf /tmp/.cache /tmp/uv-*.lock \
    && mkdir -p /tmp/.cache \
    && chmod 1777 /tmp/.cache

ENV PATH="/cli/.venv/bin:${PATH}"

# Not root. Singularity and Apptainer run the image as the caller whatever
# this says, as does Docker given `-u $(id -u):$(id -g)`. Numeric, so that a
# runtime enforcing non-root can check it.
RUN useradd --uid 10001 --user-group --create-home --no-log-init droplets
USER 10001:10001
WORKDIR /home/droplets

# No ENTRYPOINT: Nextflow requires /bin/bash to be the container entrypoint,
# so invocations spell out the command: `docker run IMAGE sc-droplets methods`.
CMD ["sc-droplets", "--help"]
