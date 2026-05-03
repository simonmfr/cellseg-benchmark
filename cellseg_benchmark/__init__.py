from importlib import import_module

from ._constants import BASE_PATH

_LAZY_SUBMODULES = {
    "adata_utils",
    "cell_annotation_utils",
    "dea_utils",
    "ficture_utils",
    "sdata_utils",
    "metrics",
}

__all__ = ["BASE_PATH", *_LAZY_SUBMODULES]


def __getattr__(name):
    if name in _LAZY_SUBMODULES:
        module = import_module(f".{name}", __name__)
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(list(globals().keys()) + list(_LAZY_SUBMODULES))