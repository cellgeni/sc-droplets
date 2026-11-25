"""Regenerate the table of methods on docs/methods/index.md from the registry.

    uv run python docs/_scripts/methods_table.py           # rewrite the table
    uv run python docs/_scripts/methods_table.py --check   # exit 1 if it is stale

tests/test_docs_are_accurate.py runs `render` and compares.
"""

from __future__ import annotations

import sys
from pathlib import Path

PAGE = Path(__file__).resolve().parent.parent / "methods" / "index.md"
START = "<!-- methods-table:start -->"
END = "<!-- methods-table:end -->"


def render() -> str:
    from sc_droplets import registry

    lines = ["| Method | Tool | Category | Versions | Page |", "|---|---|---|---|---|"]
    for cls in registry.all_methods():
        info = cls.info
        if info.status != "implemented":
            continue
        versions = "; ".join(f"{v.name}: {v.versions}" for v in info.variants) if info.versioned else "—"
        link = Path(info.docs).relative_to("methods")
        lines.append(f"| `{info.id}` | {info.tool} | {info.category} | {versions} | [{info.title}]({link}) |")
    return "\n".join(lines)


def current(text: str) -> str:
    """The table now on the page."""
    head, _, rest = text.partition(START)
    body, _, _ = rest.partition(END)
    return body.strip()


def rewrite(text: str) -> str:
    head, _, rest = text.partition(START)
    _, _, tail = rest.partition(END)
    return f"{head}{START}\n{render()}\n{END}{tail}"


def main(argv: list) -> int:
    text = PAGE.read_text()
    if START not in text or END not in text:
        print(f"{PAGE} lacks the {START} / {END} markers", file=sys.stderr)
        return 1
    if "--check" in argv:
        stale = current(text) != render()
        if stale:
            print(f"{PAGE} is stale; run docs/_scripts/methods_table.py", file=sys.stderr)
        return int(stale)
    PAGE.write_text(rewrite(text))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
