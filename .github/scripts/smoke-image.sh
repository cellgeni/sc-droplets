#!/usr/bin/env bash
# Check that an sc-droplets image runs, before anything pushes it.
#
#     bash .github/scripts/smoke-image.sh IMAGE VERSION
#
# Its help and version, that it runs unprivileged, and a call and a plot as
# the caller on a read-only root filesystem, as Singularity and Nextflow run
# it. tests.yml runs this on every pull request and before a release reaches
# PyPI; quay-on-tag.yml runs it again on the image it pushes.
set -euo pipefail

image=$1
version=$2

docker run --rm "$image" sc-droplets --help > /dev/null
reported=$(docker run --rm "$image" sc-droplets --version)
if [ "$reported" != "$version" ]; then
  echo "::error::the image reports $reported, the release is $version"
  exit 1
fi
if [ "$(docker run --rm "$image" id -u)" = 0 ]; then
  echo "::error::the image runs as root"
  exit 1
fi

# Outside the build context, so that a later build of it is not changed.
data=$(mktemp -d "${RUNNER_TEMP:-/tmp}/smoke.XXXXXX")
chmod 777 "$data"
as_caller() {
  docker run --rm --user "$(id -u):$(id -g)" --read-only --tmpfs /tmp \
    -v "$data":/data "$image" "$@"
}
as_caller python -c "import sc_droplets as scd; scd.simulate_droplets(n_cells=100, n_empty=500, n_low=1000, n_genes=100).write_h5ad('/data/sim.h5ad')"
as_caller sc-droplets inspect /data/sim.h5ad
as_caller sc-droplets call threshold /data/sim.h5ad --min-counts 150 --barcodes /data/cells.txt
as_caller sc-droplets plot knee /data/sim.h5ad -m top-n:n_cells=100 -o /data/knee.png
test -s "$data/cells.txt"
test -s "$data/knee.png"
