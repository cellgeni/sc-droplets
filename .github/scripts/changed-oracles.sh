#!/usr/bin/env bash
# Which oracle jobs a change needs, for .github/workflows/oracles.yml. Writes `tools=` with the job
# names, space-padded (` dropletutils umitools `), to $GITHUB_OUTPUT. Every job when the run was
# started by hand, when the base commit is unknown (a new branch), or when the workflow or this
# script changed; otherwise the tools whose oracles/<tool>/ changed, and those that use a changed
# shared helper.
#
#   EVENT=pull_request BASE=<sha> HEAD=<sha> bash .github/scripts/changed-oracles.sh
set -euo pipefail

all="starsolo dropletutils umitools kneedle cellsweep cellbender bustools alevinfry"
tools=""
if [ "$EVENT" = workflow_dispatch ] || [ -z "${BASE:-}" ] || ! git cat-file -e "${BASE}^{commit}" 2>/dev/null; then
  tools=$all
else
  while read -r path; do
    case "$path" in
      .github/workflows/oracles.yml | .github/scripts/changed-oracles.sh) tools="$all" ;;
      oracles/make_counts.py) tools="$tools umitools kneedle cellsweep alevinfry" ;;
      oracles/starsolo/make_mtx.py) tools="$tools starsolo dropletutils" ;;
      oracles/*/*) tools="$tools $(cut -d/ -f2 <<< "$path")" ;;
    esac
  done < <(git diff --name-only "${BASE}...${HEAD}")
fi
unique=$(tr ' ' '\n' <<< "$tools" | sed '/^$/d' | sort -u | tr '\n' ' ')
echo "tools= ${unique}" | tee -a "$GITHUB_OUTPUT"
