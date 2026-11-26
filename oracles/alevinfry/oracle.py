"""Run one alevin-fry release's program (built by harness.py) on a dataset's counts.

    python oracle.py HARNESS RELEASE COUNTS.tsv OUT_DIR EXPECTED_CELLS

COUNTS.tsv is `barcode<TAB>count` for every barcode of the matrix (oracles/make_counts.py --all);
barcode i of the program is line i. Writes OUT_DIR/cases.tsv:
`case<TAB>cells<TAB>smallest count kept<TAB>threshold` (the threshold where the rule prints one,
else -), or `case<TAB>error<TAB>message`. Every rule keeps the barcodes with count >= a threshold,
so the number of cells and the smallest count kept fix the cells.
"""

import os
import subprocess
import sys


def cases(release, expected):
    out = [("knee", ["knee"]), ("expect", ["expect", str(expected)]), ("force", ["force", str(expected)])]
    if release >= (0, 18, 0):
        out.append(("unfiltered", ["unfiltered", "10"]))
    if (0, 12, 0) <= release < (0, 18, 0):
        out += [("flex_" + name, args + ["--flex"]) for name, args in out[:3]]
    return out


def main(harness, release, counts_path, out, expected):
    version = tuple(int(p) for p in release.split("."))
    with open(counts_path) as handle:
        counts = [int(line.rstrip("\n").split("\t")[1]) for line in handle]
    text = "\n".join(str(c) for c in counts) + "\n"
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "cases.tsv"), "w") as handle:
        for name, args in cases(version, expected):
            run = subprocess.run([harness] + args, input=text, capture_output=True, text=True)
            if run.returncode != 0:
                message = next((l for l in run.stderr.splitlines() if "panicked" in l or "error" in l.lower()), run.stderr.strip()[:200])
                handle.write("%s\terror\t%s\n" % (name, message))
                print("%s %s: error %s" % (release, name, message))
                continue
            kept, threshold = [], "-"
            for line in run.stdout.splitlines():
                if line.startswith("threshold "):
                    threshold = line.split()[1]
                elif line.strip():
                    kept.append(int(line))
            smallest = min(counts[i] for i in kept) if kept else "-"
            handle.write("%s\t%d\t%s\t%s\n" % (name, len(kept), smallest, threshold))
            print("%s %s: %d cells, smallest %s, threshold %s" % (release, name, len(kept), smallest, threshold))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]))
