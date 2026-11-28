"""Check that the documentation describes the package that actually exists.

Every `sc-droplets ...` line in the README and docs is parsed against the real
command tree; every `scd.` name in a Python block must exist; every method has
its page, with the sections every method page has and every option it takes;
the generated table of methods is current; pages have the front matter the
sidebar needs; and relative links point at pages and headings that exist.
Prose drifts from code silently, and nothing else in the suite reads the docs.
"""

from __future__ import annotations

import importlib.util
import os
import re
import shlex
from pathlib import Path
from typing import Any, Dict, Iterator, List, Tuple

import pytest
from typer.main import get_command

import sc_droplets
from sc_droplets import registry
from sc_droplets.cli import app

ROOT = get_command(app)
REPO = Path(__file__).resolve().parent.parent
SITE = REPO / "docs"
PAGES = sorted(p for p in SITE.rglob("*.md") if not any(part.startswith("_") for part in p.relative_to(SITE).parts))
DOCS = [REPO / "README.md"] + PAGES

METHOD_SECTIONS = [
    "At a glance",
    "How it works",
    "Versions",
    "Parameters",
    "Run it",
    "Differences from the original",
    "Validation",
    "References",
]


def _label(path: Path) -> str:
    return str(path.relative_to(REPO))


def front_matter(text: str) -> Tuple[Dict[str, Any], str]:
    """Flat `key: value` front matter (inline `[a, b]` lists as lists) and the body."""
    if not text.startswith("---\n"):
        return {}, text
    head, _, body = text[4:].partition("\n---\n")
    found: Dict[str, Any] = {}
    for line in head.splitlines():
        pair = re.match(r"^([A-Za-z_]+):\s*(.*)$", line)
        if not pair:
            continue
        key, value = pair.group(1), pair.group(2).strip().strip('"')
        if value.startswith("[") and value.endswith("]"):
            found[key] = [v.strip() for v in value[1:-1].split(",") if v.strip()]
        elif value:
            found[key] = value
    return found, body


def _fenced(text: str, languages: str) -> Iterator[str]:
    for match in re.finditer(rf"```(?:{languages})\n(.*?)```", text, re.S):
        yield match.group(1)


# ---------------------------------------------------------------------------
# Command lines
# ---------------------------------------------------------------------------


def _invocations(path: Path) -> List[Tuple[str, str]]:
    found = []
    for block in _fenced(path.read_text(), "bash|shell|sh|console"):
        for raw in block.splitlines():
            line = re.split(r"\s+#\s", raw.strip())[0].strip()
            line = re.split(r"\s*[|>]\s*", line)[0].strip()
            if line.startswith("sc-droplets ") or line == "sc-droplets":
                found.append((line, _label(path)))
    return found


INVOCATIONS = [inv for doc in DOCS for inv in _invocations(doc)]


def test_the_docs_contain_command_lines_to_check():
    """A parsing regression here would make every case below vacuous."""
    assert len(INVOCATIONS) > 20, INVOCATIONS


def _resolve(tokens: List[str]):
    command = ROOT
    path: List[str] = []
    rest = tokens[1:]
    while rest and hasattr(command, "commands") and not rest[0].startswith("-"):
        sub = command.commands.get(rest[0])
        if sub is None:
            raise KeyError(f"{rest[0]!r} is not a command of {' '.join(['sc-droplets'] + path)}")
        path.append(rest[0])
        command = sub
        rest = rest[1:]
    return command, rest


@pytest.mark.parametrize("line, source", INVOCATIONS, ids=[f"{s}: {l[:60]}" for l, s in INVOCATIONS])
def test_every_documented_command_line_uses_real_commands_and_options(line, source):
    tokens = shlex.split(line)
    command, rest = _resolve(tokens)
    accepted = set()
    for param in command.params:
        accepted.update(param.opts)
        accepted.update(param.secondary_opts)
    accepted.add("--help")
    used = [t.split("=", 1)[0] for t in rest if t.startswith("-") and not re.fullmatch(r"-?\d+(\.\d+)?", t)]
    unknown = [t for t in used if t not in accepted]
    assert not unknown, f"{source}: {line!r} uses {unknown}, which the command does not take"


def test_every_command_has_a_section_in_the_command_reference():
    text = (SITE / "reference" / "commands.md").read_text()
    for name, command in ROOT.commands.items():
        if name in ("call", "plot"):
            for sub in command.commands:
                if name == "plot":
                    assert f"## `plot {sub}`" in text, f"no section for plot {sub}"
            continue
        assert f"## `{name}`" in text, f"no section for {name}"


# ---------------------------------------------------------------------------
# Python blocks
# ---------------------------------------------------------------------------


def _python_names() -> List[Tuple[str, str]]:
    found = []
    for path in DOCS:
        for block in _fenced(path.read_text(), "python"):
            for name in re.findall(r"\bscd\.([A-Za-z_][A-Za-z_0-9.]*)", block):
                found.append((name.rstrip("."), _label(path)))
    return found


PYTHON_NAMES = _python_names()


def test_the_docs_contain_python_to_check():
    assert len(PYTHON_NAMES) > 10


@pytest.mark.parametrize("name, source", PYTHON_NAMES, ids=[f"{s}: scd.{n}" for n, s in PYTHON_NAMES])
def test_every_documented_python_name_exists(name, source):
    target: Any = sc_droplets
    for part in name.split("."):
        assert hasattr(target, part), f"{source}: scd.{name} does not exist"
        target = getattr(target, part)


def test_the_api_reference_signatures_match_the_functions():
    text = (SITE / "reference" / "python-api.md").read_text()
    for heading in re.findall(r"^### `scd\.([\w.]+)\((.*)\)`$", text, re.M):
        name, documented = heading
        target: Any = sc_droplets
        for part in name.split("."):
            target = getattr(target, part)
        import inspect

        actual = [p for p in inspect.signature(target).parameters]
        written = [re.split(r"[=\s]", a.strip().lstrip("*"))[0] for a in documented.split(",") if a.strip() not in ("*", "...", "")]
        assert [w for w in written if w] == [a for a in actual if a in written], f"scd.{name}: documented {written}, actual {actual}"
        missing = [w for w in written if w not in actual]
        assert not missing, f"scd.{name} documents {missing}, which it does not take"


# ---------------------------------------------------------------------------
# Method pages
# ---------------------------------------------------------------------------


def _method_pages() -> Dict[str, Tuple[Path, Dict[str, Any], str]]:
    pages = {}
    for path in PAGES:
        meta, body = front_matter(path.read_text())
        for method_id in meta.get("method_ids", []):
            assert method_id not in pages, f"{method_id} is on two pages"
            pages[method_id] = (path, meta, body)
    return pages


METHOD_PAGES = _method_pages()
METHODS = [cls for cls in registry.all_methods()]


@pytest.mark.parametrize("cls", METHODS, ids=[c.info.id for c in METHODS])
def test_every_method_has_its_page_where_its_info_says(cls):
    info = cls.info
    assert info.id in METHOD_PAGES, f"no page lists method_ids: [{info.id}]"
    path, meta, _ = METHOD_PAGES[info.id]
    assert path == SITE / info.docs, f"{info.id}: MethodInfo.docs is {info.docs}, the page is {_label(path)}"
    assert meta.get("status") == info.status


@pytest.mark.parametrize("cls", METHODS, ids=[c.info.id for c in METHODS])
def test_every_method_page_has_the_sections_in_order_and_every_option(cls):
    path, _, body = METHOD_PAGES[cls.info.id]
    headings = re.findall(r"^## (.+)$", body, re.M)
    assert headings == METHOD_SECTIONS, f"{_label(path)} has sections {headings}"
    parameters = body.split("## Parameters", 1)[1].split("\n## ", 1)[0]
    for spec in cls.param_specs():
        assert f"`{spec.flag}" in parameters, f"{_label(path)} does not list {spec.flag}"


def test_every_implemented_method_page_is_for_a_registered_method():
    registered = {c.info.id for c in METHODS}
    for method_id, (path, meta, _) in METHOD_PAGES.items():
        if meta.get("status") == "implemented":
            assert method_id in registered, f"{_label(path)}: {method_id} is not registered"


#: Pages of methods not (yet) run by the package: what the method does, not how to run it here.
OTHER_SECTIONS = ["At a glance", "How it works", "Versions", "Parameters", "Status", "References"]
OTHER_PAGES = [(i, p, m) for i, (p, m, _) in METHOD_PAGES.items() if m.get("status") in ("planned", "docs-only")]
OTHER_INDEX = {"planned": SITE / "methods" / "planned" / "index.md", "docs-only": SITE / "methods" / "documented" / "index.md"}


@pytest.mark.parametrize("method_id, path, meta", OTHER_PAGES, ids=[i for i, _, _ in OTHER_PAGES])
def test_planned_and_documented_pages_have_their_sections_and_are_listed(method_id, path, meta):
    """A planned or documented-only page follows its own template, names no registered method,
    and is linked from the index of its status, wherever its tool's pages are."""
    _, body = front_matter(path.read_text())
    assert re.findall(r"^## (.+)$", body, re.M) == OTHER_SECTIONS, f"{_label(path)} has other sections"
    assert method_id not in {c.info.id for c in METHODS}, f"{_label(path)}: {method_id} is registered"
    index = OTHER_INDEX[meta["status"]]
    target = os.path.relpath(path, index.parent).replace(os.sep, "/")
    assert f"]({target})" in index.read_text(), f"{_label(index)} does not link {target}"


def test_there_are_planned_and_documented_pages():
    assert {m.get("status") for _, _, m in OTHER_PAGES} == {"planned", "docs-only"}


def test_the_methods_table_is_current():
    spec = importlib.util.spec_from_file_location("methods_table", SITE / "_scripts" / "methods_table.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    page = (SITE / "methods" / "index.md").read_text()
    assert module.current(page) == module.render(), "run: uv run python docs/_scripts/methods_table.py"


# ---------------------------------------------------------------------------
# Site structure
# ---------------------------------------------------------------------------

TITLES = {front_matter(p.read_text())[0].get("title"): p for p in PAGES}


@pytest.mark.parametrize("path", PAGES, ids=[_label(p) for p in PAGES])
def test_every_page_has_a_title_an_order_and_a_parent_that_exists(path):
    meta, _ = front_matter(path.read_text())
    assert meta.get("title"), f"{_label(path)} has no title"
    assert meta.get("nav_order"), f"{_label(path)} has no nav_order"
    for key in ("parent", "grand_parent"):
        if key in meta:
            assert meta[key] in TITLES, f"{_label(path)}: {key} {meta[key]!r} is no page's title"
    if "parent" in meta:
        parent_meta, _ = front_matter(TITLES[meta["parent"]].read_text())
        assert parent_meta.get("has_children") == "true", f"{meta['parent']} lacks has_children"


def _anchors(path: Path) -> set:
    """Heading ids as kramdown makes them."""
    _, body = front_matter(path.read_text())
    anchors = set()
    for heading in re.findall(r"^#+ (.+)$", body, re.M):
        text = re.sub(r"[`*]", "", heading).strip().lower()
        text = re.sub(r"[^\w\- ]", "", text)
        anchors.add(text.replace(" ", "-"))
    return anchors


LINKS = [
    (path, target)
    for path in DOCS
    for target in re.findall(r"\]\(([^)\s]+)\)", path.read_text())
    if not re.match(r"^(https?:|mailto:|#)", target)
]


@pytest.mark.parametrize("path, target", LINKS, ids=[f"{_label(p)} -> {t}" for p, t in LINKS])
def test_every_relative_link_points_at_a_page_and_heading_that_exist(path, target):
    file, _, anchor = target.partition("#")
    destination = (path.parent / file).resolve()
    assert destination.exists(), f"{_label(path)} links to {target}, which does not exist"
    if path.parent.is_relative_to(SITE) and destination.suffix == ".md":
        assert destination.is_relative_to(SITE), f"{_label(path)} links outside the site"
    if anchor and destination.suffix == ".md":
        assert anchor in _anchors(destination), f"{_label(path)}: no heading #{anchor} in {file}"


ASSETS = [
    (path, target)
    for path in PAGES
    for target in re.findall(
        r"\{\{\s*['\"](/[^'\"]+)['\"]\s*\|\s*relative_url\s*\}\}",
        re.sub(r"\{% raw %\}.*?\{% endraw %\}", "", path.read_text(), flags=re.S),  # examples, not links
    )
]


@pytest.mark.parametrize("path, target", ASSETS, ids=[f"{_label(p)} -> {t}" for p, t in ASSETS])
def test_every_site_asset_exists(path, target):
    """Images are linked as `{{ '/assets/...' | relative_url }}`: a relative path would break under
    the site's pretty permalinks, and the link check above skips them."""
    assert (SITE / target.lstrip("/")).is_file(), f"{_label(path)} shows {target}, which does not exist"


@pytest.mark.parametrize("path", PAGES, ids=[_label(p) for p in PAGES])
def test_inline_formulas_use_kramdown_dollars(path):
    """kramdown reads `\\(` as an escaped parenthesis, so `\\(x\\)` showed as "(x)" on the
    site; inline formulas are written `$$x$$`, which kramdown turns into KaTeX's form."""
    _, body = front_matter(path.read_text())
    prose = re.sub(r"```.*?```", "", body, flags=re.S)
    assert "\\(" not in prose, f"{_label(path)} has \\( ... \\); write $$...$$ instead"
