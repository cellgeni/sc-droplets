"""The method registry: what each method is, which versions of its tool it covers.

A method is a class derived from `Method`, decorated with `@register`. It
declares a `MethodInfo` (identity, provenance, versions) and a `Params`
dataclass (its parameters, built with `param`). The CLI builds one
`sc-droplets call <id>` subcommand per registered method from these, and the
docs test checks every method has its page.

Versions are those of the *original tool*: `cellranger@9.0` is what
Cell Ranger 9.0 did. Each `Variant` covers a half-open range of tool releases
[since, until). A `Gap` is a range whose behaviour is unknown (no public
source); asking for it is an error that names the docs page.
"""

from __future__ import annotations

import dataclasses
import enum
import functools
import re
from dataclasses import MISSING, dataclass, field
from typing import Any, ClassVar, Dict, List, Mapping, Optional, Sequence, Tuple, Type

DOCS_URL = "https://cellgeni.github.io/sc-droplets"

#: Categories, in the order the docs present them.
CATEGORIES = (
    "rank-threshold",  # a cut on the ranked totals (top-N, fixed, OrdMag)
    "knee",  # a geometric knee or inflection of the rank curve
    "density",  # a valley or cut in the distribution of log totals
    "mixture",  # a fitted mixture of distributions
    "statistical-test",  # a test of each droplet against the ambient profile
    "classifier",  # a model trained on droplet features
)

#: Count units a method works on.
UNITS = ("umis", "reads", "either")

#: Names the CLI uses for its own options; a method parameter may not take one.
RESERVED_PARAMS = frozenset(
    {
        "input",
        "output",
        "cells_only",
        "barcodes",
        "table",
        "json",
        "json_out",
        "tool_version",
        "feature_type",
        "feature_types",
        "layer",
        "genome",
        "reads",
        "seed",
        "random_state",
        "key",
        "key_added",
        "solo_feature",
        "splicing",
        "inplace",
    }
)


class UnsupportedVersion(ValueError):
    """A tool version the method cannot reproduce (a gap, or before it existed)."""


class UnknownMethod(KeyError):
    """No method is registered under that id."""


# ---------------------------------------------------------------------------
# Tool versions
# ---------------------------------------------------------------------------

_VERSION = re.compile(r"^(?:v|cellranger-)?(\d+(?:\.\d+)*)([a-z]?)$", re.IGNORECASE)


@functools.total_ordering
class ToolVersion:
    """A release of an original tool, ordered as its maintainers number them.

    Accepts `10.1.0`, `9.0`, `1.18.1`, `v0.3.2`, `cellranger-7.2.0` and STAR's
    `2.7.10a`. Missing trailing components are zeros, so `9.0 == 9.0.0`; a
    trailing letter orders after the bare number (`2.7.10 < 2.7.10a <
    2.7.10b`).
    """

    __slots__ = ("text", "_key")

    def __init__(self, text: str):
        text = str(text).strip()
        match = _VERSION.match(text)
        if not match:
            raise ValueError(
                f"{text!r} is not a tool version (expected numbers separated by dots, "
                "optionally ending in a letter, e.g. 3.0.2 or 2.7.10a)"
            )
        numbers = [int(part) for part in match.group(1).split(".")]
        while len(numbers) > 1 and numbers[-1] == 0:
            numbers.pop()
        self.text = text
        self._key = (tuple(numbers), match.group(2).lower())

    def _padded(self) -> Tuple[Tuple[int, ...], str]:
        numbers, letter = self._key
        return numbers + (0,) * (6 - len(numbers)), letter

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ToolVersion):
            return NotImplemented
        return self._padded() == other._padded()

    def __lt__(self, other: "ToolVersion") -> bool:
        return self._padded() < other._padded()

    def __hash__(self) -> int:
        return hash(self._padded())

    def __repr__(self) -> str:
        return f"ToolVersion({self.text!r})"

    def __str__(self) -> str:
        return self.text


# ---------------------------------------------------------------------------
# Method metadata
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Variant:
    """One behaviour of a method, used by tool releases [since, until).

    `settings` are flags the implementation reads; `defaults` are parameter
    defaults that differ from the `Params` dataclass for this variant.
    """

    name: str
    since: Optional[str] = None
    until: Optional[str] = None
    summary: str = ""
    settings: Mapping[str, Any] = field(default_factory=dict)
    defaults: Mapping[str, Any] = field(default_factory=dict)

    def covers(self, version: ToolVersion) -> bool:
        if self.since is not None and version < ToolVersion(self.since):
            return False
        if self.until is not None and not version < ToolVersion(self.until):
            return False
        return True

    @property
    def versions(self) -> str:
        """The range as the docs write it: `≥1.10, <1.28`, `≥2.7.8a`, `all`."""
        if self.since is None and self.until is None:
            return "all"
        if self.until is None:
            return f"≥{self.since}"
        if self.since is None:
            return f"<{self.until}"
        return f"≥{self.since}, <{self.until}"


@dataclass(frozen=True)
class Gap:
    """Tool releases [since, until) whose behaviour cannot be reproduced."""

    since: str
    until: str
    reason: str

    def covers(self, version: ToolVersion) -> bool:
        return ToolVersion(self.since) <= version < ToolVersion(self.until)


@dataclass(frozen=True)
class Reference:
    """A citation: free text plus a DOI or URL."""

    text: str
    doi: Optional[str] = None
    url: Optional[str] = None

    @property
    def link(self) -> Optional[str]:
        if self.doi:
            return f"https://doi.org/{self.doi}"
        return self.url


@dataclass(frozen=True)
class MethodInfo:
    """Everything about a method except how it computes."""

    id: str
    title: str
    tool: str
    category: str
    summary: str
    docs: str  # page path under docs/, e.g. "methods/generic/top-n.md"
    counts: str = "umis"
    deterministic: bool = True
    variants: Tuple[Variant, ...] = (Variant("default"),)
    gaps: Tuple[Gap, ...] = ()
    #: Newest release of the tool checked against; later ones get a warning.
    latest_known: Optional[str] = None
    references: Tuple[Reference, ...] = ()
    source: Optional[str] = None
    original_license: Optional[str] = None
    #: Feature types counted by default (None: every feature).
    feature_types: Optional[Tuple[str, ...]] = ("Gene Expression",)
    #: Seed the original tool uses; `random_state=None` means this one.
    default_seed: Optional[int] = None
    status: str = "implemented"

    def __post_init__(self) -> None:
        if not re.fullmatch(r"[a-z0-9]+(?:[.-][a-z0-9]+)*", self.id):
            raise ValueError(f"method id {self.id!r} must be lower-case kebab-case")
        if self.category not in CATEGORIES:
            raise ValueError(f"{self.id}: unknown category {self.category!r}")
        if self.counts not in UNITS:
            raise ValueError(f"{self.id}: unknown count unit {self.counts!r}")
        if not self.variants:
            raise ValueError(f"{self.id}: needs at least one variant")
        names = [v.name for v in self.variants]
        if len(set(names)) != len(names):
            raise ValueError(f"{self.id}: variant names repeat: {names}")

    @property
    def docs_url(self) -> str:
        page = self.docs[: -len(".md")] if self.docs.endswith(".md") else self.docs
        return f"{DOCS_URL}/{page}/"

    @property
    def versioned(self) -> bool:
        """Whether the method follows the releases of an original tool."""
        return any(v.since or v.until for v in self.variants)

    def variant(self, name: str) -> Variant:
        for v in self.variants:
            if v.name == name:
                return v
        raise KeyError(f"{self.id} has no variant {name!r}; it has {[v.name for v in self.variants]}")

    def resolve(self, tool_version: Optional[str]) -> Tuple[Variant, List[str]]:
        """The variant a tool release ran, and any warnings about the choice."""
        if tool_version is None or str(tool_version).lower() == "latest":
            return self.variants[-1], []
        if not self.versioned:
            raise UnsupportedVersion(
                f"{self.id} does not follow the releases of a tool, so it takes no "
                f"tool version (got {tool_version!r})"
            )
        version = ToolVersion(tool_version)
        for gap in self.gaps:
            if gap.covers(version):
                raise UnsupportedVersion(
                    f"{self.id}: {self.tool} {tool_version} cannot be reproduced: {gap.reason} "
                    f"See {self.docs_url}"
                )
        for v in self.variants:
            if v.covers(version):
                warnings: List[str] = []
                if self.latest_known and version > ToolVersion(self.latest_known):
                    warnings.append(
                        f"{self.tool} {tool_version} is newer than {self.latest_known}, the last "
                        f"release checked; using the {v.name} behaviour"
                    )
                return v, warnings
        first = self.variants[0]
        if first.since and version < ToolVersion(first.since):
            raise UnsupportedVersion(
                f"{self.id}: {self.tool} {tool_version} predates this method "
                f"(first in {first.since}). See {self.docs_url}"
            )
        raise UnsupportedVersion(
            f"{self.id}: no known behaviour for {self.tool} {tool_version}. See {self.docs_url}"
        )


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------


def param(
    default: Any = MISSING,
    help: str = "",
    *,
    original: Optional[str] = None,
    choices: Optional[Sequence[str]] = None,
    minimum: Optional[float] = None,
    maximum: Optional[float] = None,
) -> Any:
    """A method parameter: a dataclass field carrying its CLI help.

    `original` is the name the original tool gives it (`--expect-cells`,
    `lower`), shown on the docs page beside ours.
    """
    metadata = {
        "help": help,
        "original": original,
        "choices": tuple(choices) if choices else None,
        "minimum": minimum,
        "maximum": maximum,
    }
    if default is MISSING:
        return field(metadata=metadata)
    return field(default=default, metadata=metadata)


@dataclass(frozen=True)
class ParamSpec:
    """One parameter as the CLI and docs see it."""

    name: str
    type: type
    required: bool
    default: Any
    help: str
    original: Optional[str]
    choices: Optional[Tuple[str, ...]]
    minimum: Optional[float]
    maximum: Optional[float]

    @property
    def flag(self) -> str:
        return "--" + self.name.replace("_", "-")


def _base_type(annotation: Any) -> type:
    """int, float, str or bool from a field's annotation (which may be a string)."""
    text = annotation if isinstance(annotation, str) else getattr(annotation, "__name__", str(annotation))
    text = text.replace("typing.", "").replace(" ", "")
    for prefix, suffix in (("Optional[", "]"),):
        if text.startswith(prefix) and text.endswith(suffix):
            text = text[len(prefix) : -len(suffix)]
    text = text.replace("|None", "").replace("None|", "")
    mapping = {"int": int, "float": float, "str": str, "bool": bool}
    if text not in mapping:
        raise TypeError(f"parameter type {annotation!r} is not one of int, float, str, bool")
    return mapping[text]


class Method:
    """Base class of every method. Subclasses set `info` and `Params`, and `run`."""

    info: ClassVar[MethodInfo]
    Params: ClassVar[type]

    @classmethod
    def param_specs(cls) -> List[ParamSpec]:
        specs = []
        for f in dataclasses.fields(cls.Params):
            required = f.default is MISSING and f.default_factory is MISSING  # type: ignore[misc]
            specs.append(
                ParamSpec(
                    name=f.name,
                    type=_base_type(f.type),
                    required=required,
                    default=None if required else f.default,
                    help=f.metadata.get("help", ""),
                    original=f.metadata.get("original"),
                    choices=f.metadata.get("choices"),
                    minimum=f.metadata.get("minimum"),
                    maximum=f.metadata.get("maximum"),
                )
            )
        return specs

    @classmethod
    def build_params(cls, variant: Variant, given: Mapping[str, Any]) -> Any:
        """The Params for one run: dataclass defaults < variant defaults < given.

        `given` holds only what the caller set; None means "not set".
        """
        specs = {s.name: s for s in cls.param_specs()}
        unknown = sorted(set(given) - set(specs))
        if unknown:
            raise TypeError(
                f"{cls.info.id} has no parameter {', '.join(unknown)}; "
                f"it takes {', '.join(sorted(specs)) or 'none'}"
            )
        values: Dict[str, Any] = {}
        for name, spec in specs.items():
            value = given.get(name)
            if value is None:
                value = variant.defaults.get(name, spec.default)
            if value is None and spec.required:
                raise TypeError(f"{cls.info.id} needs {spec.flag} ({spec.help})")
            if value is not None:
                value = _coerce(cls.info.id, spec, value)
            values[name] = value
        return cls.Params(**values)

    def run(self, droplets: Any, params: Any, variant: Variant, random_state: Optional[int]) -> Any:
        raise NotImplementedError


_BOOLEAN_WORDS = {"true": True, "yes": True, "1": True, "false": False, "no": False, "0": False}


def _coerce(method: str, spec: ParamSpec, value: Any) -> Any:
    """Check and convert one parameter value, with a message naming it."""
    if isinstance(value, enum.Enum):
        value = value.value  # a choice from the CLI (a str Enum, so check this first)
    try:
        if spec.type is bool:
            if isinstance(value, str) and value.lower() in _BOOLEAN_WORDS:
                value = _BOOLEAN_WORDS[value.lower()]
            if not isinstance(value, bool):
                raise TypeError
        elif spec.type is int:
            if isinstance(value, bool) or (isinstance(value, float) and not value.is_integer()):
                raise TypeError
            value = int(value)
        elif spec.type is float:
            if isinstance(value, bool):
                raise TypeError
            value = float(value)
        else:
            value = str(value)
    except (TypeError, ValueError):
        raise TypeError(f"{method}: {spec.flag} must be {spec.type.__name__}, got {value!r}") from None
    if spec.choices and value not in spec.choices:
        raise ValueError(f"{method}: {spec.flag} must be one of {', '.join(spec.choices)}, got {value!r}")
    if spec.minimum is not None and value < spec.minimum:
        raise ValueError(f"{method}: {spec.flag} must be at least {spec.minimum}, got {value!r}")
    if spec.maximum is not None and value > spec.maximum:
        raise ValueError(f"{method}: {spec.flag} must be at most {spec.maximum}, got {value!r}")
    return value


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

_REGISTRY: Dict[str, Type[Method]] = {}


def register(cls: Type[Method]) -> Type[Method]:
    """Class decorator: add a method to the registry."""
    info = cls.info
    if info.id in _REGISTRY and _REGISTRY[info.id] is not cls:
        raise ValueError(f"method id {info.id!r} is registered twice")
    if not dataclasses.is_dataclass(cls.Params):
        raise TypeError(f"{info.id}: Params must be a dataclass")
    for spec in cls.param_specs():
        if spec.name in RESERVED_PARAMS:
            raise ValueError(f"{info.id}: parameter {spec.name!r} clashes with a CLI option")
    for variant in info.variants:
        names = {s.name for s in cls.param_specs()}
        stray = set(variant.defaults) - names
        if stray:
            raise ValueError(f"{info.id}: variant {variant.name} sets defaults for unknown {stray}")
    _REGISTRY[info.id] = cls
    return cls


def _load_methods() -> None:
    """Import every method module, so each registers itself."""
    from sc_droplets import methods  # noqa: F401  (imports its submodules)


def get(method_id: str) -> Type[Method]:
    _load_methods()
    try:
        return _REGISTRY[method_id]
    except KeyError:
        known = ", ".join(sorted(_REGISTRY))
        raise UnknownMethod(f"no method {method_id!r}; known methods: {known}") from None


def all_methods() -> List[Type[Method]]:
    """Every registered method, by tool then id."""
    _load_methods()
    return sorted(_REGISTRY.values(), key=lambda c: (c.info.tool.lower(), c.info.id))


@dataclass(frozen=True)
class Spec:
    """A method spec: `id[@tool_version][:name=value,...]`.

    `top-n:n_cells=3000`, `dropletutils-emptydrops@1.26:lower=200`. Values
    stay strings here; `Method.build_params` converts them.
    """

    id: str
    version: Optional[str] = None
    params: Mapping[str, str] = field(default_factory=dict)


def parse_spec(spec: str) -> Spec:
    """Parse a method spec `id[@tool_version][:name=value,...]`."""
    text = spec.strip()
    head, _, tail = text.partition(":")
    params: Dict[str, str] = {}
    if tail:
        for item in tail.split(","):
            name, sep, value = item.partition("=")
            if not sep or not name.strip():
                raise ValueError(f"{spec!r}: parameters after ':' are name=value, separated by commas")
            params[name.strip().replace("-", "_")] = value.strip()
    method_id, at, version = head.partition("@")
    if at and not version:
        raise ValueError(f"{spec!r}: nothing after '@'")
    if not method_id:
        raise ValueError(f"{spec!r}: no method id")
    return Spec(method_id.strip(), version.strip() or None, params)


def default_key(method_id: str, tool_version: Optional[str] = None, params: Optional[Mapping[str, Any]] = None) -> str:
    """The obs/uns key a call is stored under: `cellranger_emptydrops_9_0`.

    Parameters given in a spec are part of the key, so two calls of one
    method with different settings do not overwrite each other.
    """
    parts = [method_id]
    if tool_version:
        parts.append(str(tool_version))
    for name, value in sorted((params or {}).items()):
        parts.append(f"{name}_{value}")
    key = re.sub(r"[^a-z0-9]+", "_", "_".join(parts).lower())
    return key.strip("_")
