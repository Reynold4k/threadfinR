"""The Python package ``threadfin`` as the R package exports it (``tf``): the same objects, R-friendly argument types.

Attribute access is forwarded to :mod:`threadfin`. Submodules come back as proxies of themselves, and Threadfin's
functions are wrapped so that values from R arrive in the type the Python signature asks for: a whole number such
as ``random_state = 0`` or ``n_jobs = 8`` (an R double) becomes a Python ``int`` where the parameter takes an
integer, a vector of whole numbers becomes a sequence of ``int`` where the parameter takes one, and a single string
reaches a parameter that takes only a sequence of strings as a one-element list. Nothing else changes: each wrapper
calls the Threadfin function itself and keeps its name, docstring and signature, so results are identical to
calling Threadfin from Python. Threadfin's warnings are printed on one line each. Only underscore names are defined
here, so no Threadfin name is shadowed.
"""

import functools as _functools
import sys as _sys
import types as _types

import threadfin as _target

if getattr(_target, "__file__", None) is None:
    # A folder named "threadfin" without the package (a namespace package), e.g. R's own installed package folder
    # in the working directory: not the Python package.
    _where = ", ".join(str(p) for p in getattr(_target, "__path__", []))
    del _sys.modules["threadfin"]
    raise ModuleNotFoundError(f"No module named 'threadfin' (only a folder of that name was found: {_where})")

from threadfin_r.glue import call_recording as _call_recording
from threadfin_r.glue import coerce_call as _coerce_call

_WRAPPED = {}
_PROXIES = {}


def _is_threadfin(obj):
    return str(getattr(obj, "__module__", "") or "").split(".")[0] == "threadfin"


def _wrap(fn):
    """``fn`` with R values coerced to the types of its signature (one wrapper per function)."""
    hit = _WRAPPED.get(id(fn))
    if hit is not None and hit.__wrapped__ is fn:
        return hit

    @_functools.wraps(fn)
    def call(*args, **kwargs):
        args, kwargs = _coerce_call(fn, args, kwargs)
        return _call_recording(fn, args, kwargs, keep=False)

    _WRAPPED[id(fn)] = call
    return call


def _forward(target, name):
    value = getattr(target, name)
    if isinstance(value, _types.ModuleType):
        return _proxy(value) if value.__name__.split(".")[0] == "threadfin" else value
    if callable(value) and not isinstance(value, type) and _is_threadfin(value):
        return _wrap(value)
    return value


def _proxy(target):
    """A module whose attributes are those of ``target`` (functions wrapped, submodules proxied)."""
    hit = _PROXIES.get(target.__name__)
    if hit is not None:
        return hit
    mod = _types.ModuleType(target.__name__, target.__doc__)

    def __getattr__(name):  # PEP 562
        return _forward(target, name)

    def __dir__():
        return dir(target)

    mod.__getattr__ = __getattr__
    mod.__dir__ = __dir__
    mod.__threadfin__ = target
    _PROXIES[target.__name__] = mod
    return mod


def __getattr__(name):  # PEP 562: tf$run, tf$tl, tf$`__version__`, ...
    return _forward(_target, name)


def __dir__():
    return dir(_target)
