"""Every method, one module per tool page.

Importing this package imports each module, and each module's `@register`
adds its methods to the registry. A new tool's module goes in `_MODULES`.
"""

from importlib import import_module

_MODULES = ("generic", "cellranger", "starsolo", "dropletutils", "umitools", "kneedle")

for _name in _MODULES:
    import_module(f"{__name__}.{_name}")
