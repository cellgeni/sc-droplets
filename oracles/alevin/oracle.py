"""Run one salmon release's alevin program (built from harness.py) on a dataset's counts.

    python oracle.py PROGRAM RELEASE COUNTS.tsv OUT_DIR EXPECTED_CELLS

COUNTS.tsv is `barcode<TAB>count` for the barcodes with a count, in the matrix's order
(oracles/make_counts.py); barcode i of the program is line i. Writes OUT_DIR/cases.tsv:
`case<TAB>name<TAB>value` (knee, gauss, gauss_taken, cells, low, whitelist, cells_min,
cells_above, whitelist_digest), or `case<TAB>error<TAB>message`.

alevin's high-confidence barcodes are its first `cells` in std::sort's order; where counts tie
at the boundary, which are taken depends on that order. So the program's cells are recorded as
the smallest count among them (cells_min) and how many have more (cells_above): with the count of
cells, these fix the cells up to the tied barcodes. The whitelist is recorded whole, as a count
and the SHA-256 of its sorted barcodes.
"""

import hashlib
import os
import subprocess
import sys


def cases(release, expected):
    version = tuple(int(p) for p in release.split("."))
    out = [("knee", ["knee"])]
    if version >= (0, 12, 0):
        out += [("expect", ["expect", str(expected)]), ("force", ["force", str(expected)])]
    if version >= (0, 13, 0):
        out.append(("keep", ["keep", "0.005"]))
    return out


def main(program, release, counts_path, out, expected):
    barcodes, counts = [], []
    with open(counts_path) as handle:
        for line in handle:
            barcode, count = line.rstrip("\n").split("\t")
            barcodes.append(barcode)
            counts.append(int(count))
    text = "\n".join(str(c) for c in counts) + "\n"
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "cases.tsv"), "w") as handle:
        for name, args in cases(release, expected):
            run = subprocess.run([program] + args, input=text, capture_output=True, text=True)
            if run.returncode != 0:
                message = (run.stderr.strip() or run.stdout.strip().splitlines()[-1:] or ["exit %d" % run.returncode])
                message = message if isinstance(message, str) else message[0]
                handle.write("%s\terror\texit %d: %s\n" % (name, run.returncode, message.splitlines()[0][:200]))
                print("%s %s: exit %d" % (release, name, run.returncode))
                continue
            values, cells, white = {}, [], []
            for line in run.stdout.splitlines():
                key, *rest = line.split()
                if key == "cell":
                    cells.append(int(rest[0]))
                elif key == "white":
                    white.append(int(rest[0]))
                elif key == "gauss":
                    values["gauss"], values["gauss_taken"] = rest[0], rest[1]
                else:
                    values[key] = rest[0]
            smallest = min(counts[i] for i in cells) if cells else 0
            values["cells_min"] = smallest
            values["cells_above"] = sum(1 for i in cells if counts[i] > smallest)
            names = sorted(barcodes[i] for i in white)
            values["whitelist"] = len(names)
            values["whitelist_digest"] = hashlib.sha256("\n".join(names).encode()).hexdigest()
            for key, value in values.items():
                handle.write("%s\t%s\t%s\n" % (name, key, value))
            print("%s %s: %s" % (release, name, {k: v for k, v in values.items() if k != "whitelist_digest"}))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]))
