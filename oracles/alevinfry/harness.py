"""Generate a small Rust program from one alevin-fry release's own source, for the oracle.

    python harness.py SRC_DIR OUT_DIR

SRC_DIR is a checkout of alevin-fry at a tag. The program takes alevin-fry's knee functions
(`struct Point`, `distance_to_line`, `get_max_distance_index`, `get_knee`), its
`CellFilterMethod` and its rules verbatim:

- 0.18.0 and later: `select_retained_barcodes` as it is;
- earlier: the `match filter_meth { ... }` block that turns the sorted frequencies into a
  threshold and a permit list, inside a function giving it the variables it uses;
- 0.12.0 - 0.17.1, multi-barcode input: the statements of the knee and the other rules up to
  `let threshold = ...;`.

Only names are adapted: module paths (`knee_finding::`, `libradicl::`) are dropped and hasher type
parameters of `HashMap` become the default, which changes no arithmetic. The program reads one
count per line (barcode i is line i) and prints the kept barcodes, one per line, then
`threshold <t>` where the rule computes one.

    harness MODE [VALUE] [--flex] < counts.txt
    MODE: knee | expect E | force F | unfiltered M
"""

import os
import re
import sys


def read_sources(src):
    texts = []
    for root, _, files in os.walk(src):
        for name in files:
            if name in ("cellfilter.rs", "knee_finding.rs", "lib.rs") and "atac" not in root:
                with open(os.path.join(root, name)) as handle:
                    texts.append(handle.read())
    return "\n".join(texts)


def block(text, start):
    """The text from `start` to the brace closing the first `{` after it."""
    i = text.index("{", start)
    depth = 0
    for j in range(i, len(text)):
        if text[j] == "{":
            depth += 1
        elif text[j] == "}":
            depth -= 1
            if depth == 0:
                return text[start:j + 1]
    raise ValueError("unbalanced braces")


def item(text, pattern):
    m = re.search(pattern, text)
    return block(text, m.start()) if m else None


def adapt(code):
    code = re.sub(r"HashMap<u64,\s*u64,\s*[^>]*?(<[^>]*>)?>", "HashMap<u64, u64>", code)
    code = code.replace("knee_finding::get_knee", "get_knee")
    code = code.replace("libradicl::", "")
    code = code.replace("afutils::", "")
    return code


def statements_until_threshold(arm):
    """The arm's statements up to and including `let threshold = ...;`."""
    body = arm[arm.index("{") + 1:]
    depth, at = 0, None
    for m in re.finditer(r"let threshold\s*=|[{}()\[\]]", body):
        token = m.group()
        if token in "{([":
            depth += 1
        elif token in "})]":
            depth -= 1
        elif depth == 0:  # the arm's own statement, not one nested in a match
            at = m
            break
    depth, j = 0, at.end()
    while True:
        c = body[j]
        depth += c in "({["
        depth -= c in ")}]"
        if c == ";" and depth == 0:
            return body[:j + 1]
        j += 1


def main(src, out):
    text = read_sources(src)
    point = item(text, r"\n(pub )?struct Point\b")
    knee_fns = [item(text, r"\n\s*(pub(\(crate\))? )?fn " + name + r"\b") for name in ("distance_to_line", "get_max_distance_index", "get_knee")]
    enum = item(text, r"\npub enum CellFilterMethod\b")
    enum = "#[derive(Clone, Debug)]\n" + enum
    parts = [point, *knee_fns, enum]
    select = item(text, r"\n\s*fn select_retained_barcodes\b")
    if select:
        parts.append(select)
        single = """
fn single(hist: HashMap<u64, u64>, method: CellFilterMethod, min_reads: u64, log: &slog::Logger) -> Vec<u64> {
    select_retained_barcodes(&hist, &method, min_reads, 16, log)
}
fn permit_list_from_file<P>(_ifile: P, _bclen: u16) -> Vec<u64> { unimplemented!() }"""
    else:
        start = text.index("let mut freq: Vec<u64> = hm.values()")
        head_start = [m.start() for m in re.finditer(r"\n(pub(\(crate\))? )?fn ", text[:start])][-1]
        head = text[head_start:start]
        by_ref = re.search(r"filter_meth:\s*&", head) is not None
        hm_ref = re.search(r"hm:\s*&", head) is not None
        log_ref = re.search(r"log:\s*&", head) is not None
        match = block(text, text.index("match filter_meth", start))
        threshold_fn = item(text, r"\npub fn permit_list_from_threshold\b")
        parts.append(threshold_fn)
        single = f"""
fn single(hist: HashMap<u64, u64>, method: CellFilterMethod, _min_reads: u64, logger: &slog::Logger) -> Vec<u64> {{
    let hm = {"&hist" if hm_ref else "hist.clone()"};
    let filter_meth = {"&method" if by_ref else "method.clone()"};
    let log = {"logger" if log_ref else "logger.clone()"};
    #[allow(unused_variables)]
    let ft_vals = FtVals {{ bclen: 16 }};
    #[allow(unused_variables)]
    let barcode_len: u16 = 16;
    let mut freq: Vec<u64> = hist.values().cloned().collect();
    freq.sort_unstable();
    freq.reverse();
    let valid_bc: Vec<u64>;
    {match}
    valid_bc
}}
struct FtVals {{ bclen: u16 }}
fn permit_list_from_file<P>(_ifile: P, _bclen: u16) -> Vec<u64> {{ unimplemented!() }}"""
    parts.append(single)
    multi = re.search(r"fn do_generate_permit_list_multi_bc\b", text)
    flex = """
fn flex(_cell_hist: HashMap<u64, u64>, _method: CellFilterMethod, _log: &slog::Logger) -> u64 { panic!("no multi-barcode rules in this release") }"""
    if multi and not select:
        body = block(text, multi.start())
        knee_arm = block(body, body.index("CellFilterMethod::KneeFinding =>", body.index("cell_hist")))
        other_arm = block(body, body.index("_ =>", body.index("CellFilterMethod::KneeFinding =>", body.index("cell_hist"))))
        flex = f"""
struct Opts {{ fmeth: CellFilterMethod }}
fn flex(cell_hist: HashMap<u64, u64>, method: CellFilterMethod, log: &slog::Logger) -> u64 {{
    let gpl_opts = Opts {{ fmeth: method }};
    match &gpl_opts.fmeth {{
        CellFilterMethod::KneeFinding => {{ {statements_until_threshold(knee_arm)} threshold }}
        _ => {{ {statements_until_threshold(other_arm)} threshold }}
    }}
}}"""
    parts.append(flex)
    main_rs = """#![allow(dead_code, unused_imports, unused_mut, clippy::all)]
use slog::{info, warn, Drain};
use std::cmp;
use std::collections::HashMap;
use std::io::{self, BufRead};
use std::path::PathBuf;

""" + adapt("\n\n".join(parts)) + """

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let log = slog::Logger::root(slog::Discard, slog::o!());
    let mut hist: HashMap<u64, u64> = HashMap::new();
    for (i, line) in io::stdin().lock().lines().enumerate() {
        let count: u64 = line.unwrap().trim().parse().unwrap();
        if count > 0 {
            hist.insert(i as u64, count);
        }
    }
    let value: usize = args.get(2).filter(|v| !v.starts_with("--")).map(|v| v.parse().unwrap()).unwrap_or(0);
    let method = match args[1].as_str() {
        "knee" => CellFilterMethod::KneeFinding,
        "expect" => CellFilterMethod::ExpectCells(value),
        "force" => CellFilterMethod::ForceCells(value),
        UNFILTERED
        other => panic!("unknown mode {}", other),
    };
    if args.iter().any(|a| a == "--flex") {
        let threshold = flex(hist.clone(), method, &log);
        let mut kept: Vec<u64> = hist.iter().filter(|(_, &c)| c >= threshold).map(|(&b, _)| b).collect();
        kept.sort_unstable();
        for b in kept { println!("{}", b); }
        println!("threshold {}", threshold);
    } else {
        let mut kept = single(hist, method, value as u64, &log);
        kept.sort_unstable();
        for b in kept { println!("{}", b); }
    }
}
"""
    unfiltered = ('"unfiltered" => CellFilterMethod::UnfilteredExternalList(PathBuf::from("-"), value),'
                  if "UnfilteredExternalList" in enum else "")
    main_rs = main_rs.replace("UNFILTERED", unfiltered)
    os.makedirs(os.path.join(out, "src"), exist_ok=True)
    with open(os.path.join(out, "src", "main.rs"), "w") as handle:
        handle.write(main_rs)
    with open(os.path.join(out, "Cargo.toml"), "w") as handle:
        handle.write('[package]\nname = "harness"\nversion = "0.1.0"\nedition = "2021"\n\n'
                     '[dependencies]\nslog = "2"\n\n[profile.release]\nopt-level = 3\n')


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
