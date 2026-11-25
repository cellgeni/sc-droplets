"""Tool versions, variant resolution, parameters and specs."""

from __future__ import annotations

from dataclasses import dataclass

import pytest

from sc_droplets import registry
from sc_droplets.registry import (
    Gap,
    Method,
    MethodInfo,
    ToolVersion,
    UnsupportedVersion,
    Variant,
    param,
)


@pytest.mark.parametrize(
    "low, high",
    [
        ("2.7.9a", "2.7.10a"),
        ("2.7.10", "2.7.10a"),
        ("2.7.10a", "2.7.10b"),
        ("1.18.0", "1.18.1"),
        ("v0.3.2", "0.4.0"),
        ("cellranger-7.2.0", "8.0"),
        ("9.0.1", "10.0"),
    ],
)
def test_tool_versions_order_as_their_maintainers_number_them(low, high):
    assert ToolVersion(low) < ToolVersion(high)
    assert not ToolVersion(high) < ToolVersion(low)


def test_trailing_zeros_and_prefixes_do_not_change_a_version():
    assert ToolVersion("9.0") == ToolVersion("9.0.0") == ToolVersion("cellranger-9.0.0")
    assert ToolVersion("v0.3.2") == ToolVersion("0.3.2")
    assert hash(ToolVersion("9.0")) == hash(ToolVersion("9"))


@pytest.mark.parametrize("bad", ["", "latest", "1.x", "2.7.10ab", "-1"])
def test_a_malformed_version_is_refused_with_an_example(bad):
    with pytest.raises(ValueError, match="e.g. 3.0.2"):
        ToolVersion(bad)


def _info(**kwargs) -> MethodInfo:
    defaults = dict(
        id="demo-method",
        title="Demo",
        tool="Demo",
        category="knee",
        summary="A method for the tests.",
        docs="methods/demo.md",
        variants=(
            Variant("v1", since="1.0", until="2.0", summary="first"),
            Variant("v2", since="2.0", until="4.0", summary="second"),
            Variant("v4", since="5.0", summary="third"),
        ),
        gaps=(Gap("4.0", "5.0", "no public source for these releases."),),
        latest_known="5.2",
    )
    defaults.update(kwargs)
    return MethodInfo(**defaults)


def test_each_release_resolves_to_the_variant_that_covers_it():
    info = _info()
    assert info.resolve("1.0")[0].name == "v1"
    assert info.resolve("1.9.9")[0].name == "v1"
    assert info.resolve("2.0")[0].name == "v2"
    assert info.resolve("5.1")[0].name == "v4"


def test_no_version_means_the_latest_behaviour():
    info = _info()
    assert info.resolve(None) == (info.variants[-1], [])
    assert info.resolve("latest")[0].name == "v4"


def test_a_release_in_a_gap_is_an_error_naming_the_page():
    with pytest.raises(UnsupportedVersion, match="cannot be reproduced.*methods/demo/"):
        _info().resolve("4.5")


def test_a_release_before_the_method_existed_is_an_error():
    with pytest.raises(UnsupportedVersion, match="predates"):
        _info().resolve("0.9")


def test_a_release_newer_than_any_checked_warns_and_uses_the_latest():
    variant, warnings = _info().resolve("6.0")
    assert variant.name == "v4"
    assert warnings and "newer than 5.2" in warnings[0]


def test_a_method_without_versions_refuses_a_tool_version():
    info = _info(variants=(Variant("default"),), gaps=(), latest_known=None)
    assert not info.versioned
    with pytest.raises(UnsupportedVersion, match="takes no tool version"):
        info.resolve("1.0")


def test_variant_ranges_are_written_as_the_docs_write_them():
    info = _info()
    assert [v.versions for v in info.variants] == ["≥1.0, <2.0", "≥2.0, <4.0", "≥5.0"]
    assert Variant("x").versions == "all"
    assert Variant("x", until="3").versions == "<3"


@pytest.mark.parametrize(
    "kwargs, message",
    [
        (dict(id="Bad_Id"), "kebab-case"),
        (dict(category="vibes"), "unknown category"),
        (dict(counts="molecules"), "unknown count unit"),
        (dict(variants=()), "at least one variant"),
        (dict(variants=(Variant("a"), Variant("a"))), "repeat"),
    ],
)
def test_method_info_refuses_inconsistent_metadata(kwargs, message):
    with pytest.raises(ValueError, match=message):
        _info(**kwargs)


def test_the_docs_url_drops_the_md_suffix():
    assert _info().docs_url == "https://cellgeni.github.io/sc-droplets/methods/demo/"


class _Demo(Method):
    info = _info()

    @dataclass
    class Params:
        n: int = param(help="how many", minimum=1)
        low: float = param(15.0, "a threshold")
        mode: str = param("a", "a choice", choices=("a", "b"))
        strict: bool = param(False, "a flag")


def test_parameters_take_dataclass_then_variant_then_given_values():
    variant = Variant("v", defaults={"low": 5.0})
    built = _Demo.build_params(variant, {"n": 3})
    assert (built.n, built.low, built.mode, built.strict) == (3, 5.0, "a", False)
    built = _Demo.build_params(variant, {"n": 3, "low": 7, "mode": "b", "strict": True})
    assert (built.low, built.mode, built.strict) == (7.0, "b", True)


def test_spec_strings_are_converted_to_the_parameter_types():
    built = _Demo.build_params(Variant("v"), {"n": "4", "low": "2.5", "strict": "yes"})
    assert (built.n, built.low, built.strict) == (4, 2.5, True)


@pytest.mark.parametrize(
    "given, error, message",
    [
        ({}, TypeError, "needs --n"),
        ({"n": 2, "colour": "red"}, TypeError, "no parameter colour"),
        ({"n": 2.5}, TypeError, "--n must be int"),
        ({"n": 0}, ValueError, "at least 1"),
        ({"n": 2, "mode": "c"}, ValueError, "one of a, b"),
        ({"n": 2, "strict": "maybe"}, TypeError, "--strict must be bool"),
    ],
)
def test_bad_parameters_are_refused_naming_the_option(given, error, message):
    with pytest.raises(error, match=message):
        _Demo.build_params(Variant("v"), given)


def test_param_specs_describe_each_parameter():
    specs = {s.name: s for s in _Demo.param_specs()}
    assert specs["n"].required and specs["n"].type is int and specs["n"].flag == "--n"
    assert not specs["low"].required and specs["low"].default == 15.0
    assert specs["mode"].choices == ("a", "b")


def test_registering_a_parameter_named_like_a_cli_option_is_refused():
    class Clash(Method):
        info = _info(id="clash-method")

        @dataclass
        class Params:
            seed: int = param(1, "clashes with --seed")

    with pytest.raises(ValueError, match="clashes with a CLI option"):
        registry.register(Clash)


def test_a_variant_default_for_an_unknown_parameter_is_refused():
    class Stray(Method):
        info = _info(id="stray-method", variants=(Variant("v", defaults={"nope": 1}),), gaps=())

        @dataclass
        class Params:
            n: int = param(1, "n")

    with pytest.raises(ValueError, match="unknown"):
        registry.register(Stray)


def test_an_unknown_method_lists_the_known_ones():
    with pytest.raises(registry.UnknownMethod, match="top-n"):
        registry.get("no-such-method")


@pytest.mark.parametrize(
    "text, expected",
    [
        ("top-n", ("top-n", None, {})),
        ("dropletutils-emptydrops@1.26", ("dropletutils-emptydrops", "1.26", {})),
        ("top-n:n_cells=3000", ("top-n", None, {"n_cells": "3000"})),
        ("x@2.7.10a:a=1, b-c=two", ("x", "2.7.10a", {"a": "1", "b_c": "two"})),
    ],
)
def test_specs_split_into_id_version_and_parameters(text, expected):
    spec = registry.parse_spec(text)
    assert (spec.id, spec.version, dict(spec.params)) == expected


@pytest.mark.parametrize("bad", ["top-n@", "@1.0", "top-n:n_cells", "top-n:=3"])
def test_malformed_specs_are_refused(bad):
    with pytest.raises(ValueError):
        registry.parse_spec(bad)


def test_default_keys_carry_the_version_and_spec_parameters():
    assert registry.default_key("cellranger-emptydrops") == "cellranger_emptydrops"
    assert registry.default_key("cellranger-emptydrops", "9.0") == "cellranger_emptydrops_9_0"
    assert registry.default_key("top-n", None, {"n_cells": "3000"}) == "top_n_n_cells_3000"


def test_all_methods_are_sorted_by_tool_then_id():
    ids = [(c.info.tool.lower(), c.info.id) for c in registry.all_methods()]
    assert ids == sorted(ids)
    assert {"top-n", "threshold"} <= {i for _, i in ids}
