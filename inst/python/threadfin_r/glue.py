"""Glue between R and the Python package Threadfin: argument types, objects built from R data, results and figures.

No analysis happens here. Every function either moves data between the two languages without changing a value or
draws a Threadfin figure to a file.
"""

from __future__ import annotations

import ast
import inspect
import io
import os
import sys
import warnings
from typing import Any

import numpy as np
import pandas as pd

__all__ = [
    "anndata_from_r",
    "attr",
    "call_path",
    "call_recording",
    "coerce_call",
    "contigs_from_r",
    "frame_for_r",
    "frame_from_r",
    "frame_value",
    "item",
    "keys",
    "last_warnings",
    "matrix_for_r",
    "method",
    "read_h5ad",
    "read_png",
    "render",
    "result_from_adata",
    "run_columns",
    "series_from_r",
    "value_for_r",
    "write_text",
]

_SEQUENCE = ("Sequence", "list", "tuple", "Iterable", "Collection", "set", "frozenset", "Set")
_INT = ("int", "RngLike", "np.integer", "numpy.integer")
_STRING = ("str", "Any", "Path", "os.PathLike", "PathLike")


# ------------------------------------------------------------------------------------------------- argument types


def _split(text: str, sep: str) -> list[str]:
    """Split at ``sep`` outside brackets."""
    out, depth, cur = [], 0, []
    for ch in text:
        if ch in "[(":
            depth += 1
        elif ch in "])":
            depth -= 1
        if ch == sep and depth == 0:
            out.append("".join(cur).strip())
            cur = []
        else:
            cur.append(ch)
    out.append("".join(cur).strip())
    return [m for m in out if m]


def _members(annotation: Any) -> list[str]:
    """Top-level members of a union annotation (``"int | None"`` -> ``["int", "None"]``; Optional / Union expanded)."""
    if annotation is inspect.Parameter.empty:
        return []
    text = annotation if isinstance(annotation, str) else inspect.formatannotation(annotation)
    out: list[str] = []
    for m in _split(text, "|"):
        if m.startswith(("Optional[", "typing.Optional[")) and m.endswith("]"):
            out += _members(m[m.index("[") + 1 : -1]) + ["None"]
        elif m.startswith(("Union[", "typing.Union[")) and m.endswith("]"):
            for part in _split(m[m.index("[") + 1 : -1], ","):
                out += _members(part)
        else:
            out.append(m)
    return out


def _literals(members: list[str]) -> set[Any]:
    vals: set[Any] = set()
    for m in members:
        if m.startswith(("Literal[", "typing.Literal[")) and m.endswith("]"):
            try:
                v = ast.literal_eval("(" + m[m.index("[") + 1 : -1] + ",)")
            except (ValueError, SyntaxError):
                continue
            vals.update(x for x in v if isinstance(x, str | int | float | bool))
    return vals


def _whole(value: Any) -> bool:
    return isinstance(value, float) and not isinstance(value, bool) and value.is_integer()


def _coerce(param: inspect.Parameter, value: Any) -> Any:
    members = _members(param.annotation)
    default = param.default
    if _whole(value):
        takes_int = any(m in _INT for m in members) or (
            not members and isinstance(default, int) and not isinstance(default, bool)
        )
        takes_float = "float" in members or (not members and isinstance(default, float))
        if takes_int and not takes_float:
            return int(value)
    if isinstance(value, list | tuple) and value and all(_whole(v) or isinstance(v, int) for v in value):
        # a vector of whole numbers where the parameter takes a sequence of integers, e.g. size_strata=c(2, 5, 20)
        if any(m.split("[")[0] in _SEQUENCE and "[int" in m.replace(" ", "") for m in members):
            return type(value)(int(v) for v in value)
    if isinstance(value, str):
        if members:
            sequence = any(m.split("[")[0] in _SEQUENCE for m in members)
            single = any(m.split("[")[0] in _STRING for m in members) or value in _literals(members)
            if sequence and not single:
                return [value]
        elif isinstance(default, list | tuple) and not isinstance(default, str):
            # an unannotated parameter whose default is a sequence (e.g. within=()): one name is a list of one
            return [value]
    if isinstance(value, dict) and any(m.startswith("dict[") and "list[" in m for m in members):
        # dict[str, list[str]] (gene sets): a single gene is a list of one gene
        return {k: [v] if isinstance(v, str) else v for k, v in value.items()}
    return value


def coerce_call(fn: Any, args: tuple[Any, ...], kwargs: dict[str, Any]) -> tuple[tuple[Any, ...], dict[str, Any]]:
    """Arguments from R in the types ``fn`` takes (see the module documentation of :mod:`threadfin_r`).

    Only parameters whose annotation (or default) asks for an integer, a sequence or a dict of lists are touched;
    when the signature cannot be read or the arguments do not bind, everything is passed on unchanged so that the
    function raises its own error.
    """
    try:
        sig = inspect.signature(fn)
        bound = sig.bind_partial(*args, **kwargs)
    except (TypeError, ValueError):
        return args, kwargs
    for name in list(bound.arguments):
        p = sig.parameters[name]
        if p.kind in (inspect.Parameter.VAR_POSITIONAL, inspect.Parameter.VAR_KEYWORD):
            continue
        bound.arguments[name] = _coerce(p, bound.arguments[name])
    return tuple(bound.args), dict(bound.kwargs)


# -------------------------------------------------------------------------------------------- object access


_LAST_WARNINGS: list[str] = []


def _threadfin_dirs() -> tuple[str, ...]:
    import threadfin

    here = os.path.dirname(os.path.abspath(__file__))
    return tuple(os.path.normcase(os.path.abspath(p)) for p in (os.path.dirname(threadfin.__file__), here))


def _threadfin_warning(w: Any) -> bool:
    """A warning Threadfin itself gave: a UserWarning raised in Threadfin's code or, through ``stacklevel``,
    attributed to its caller (this glue, or no Python frame at all when R called the function)."""
    if w.category is not UserWarning:
        return False
    filename = str(w.filename or "")
    if filename in ("sys", "<sys>"):              # stacklevel beyond the Python frames (R called the function)
        return True
    path = os.path.normcase(os.path.abspath(filename))
    return any(path.startswith(d + os.sep) for d in _threadfin_dirs())


def call_recording(fn: Any, args: Any, kwargs: Any, *, keep: bool = True) -> Any:
    """``fn(*args, **kwargs)`` with Threadfin's warnings kept for R (``keep``) or printed on one line each.

    Warnings of other packages are shown as Python shows them. The warning filters in force apply as usual.
    """
    if keep:
        _LAST_WARNINGS.clear()
    with warnings.catch_warnings(record=True) as caught:
        out = fn(*tuple(args or ()), **dict(kwargs or {}))
    for w in caught:
        if _threadfin_warning(w):
            text = str(w.message)
            if keep:
                _LAST_WARNINGS.append(text)
            else:
                print(f"UserWarning: {text}", file=sys.stderr)
        else:
            warnings.showwarning(w.message, w.category, w.filename, w.lineno, w.file, w.line)
    return out


def last_warnings() -> list[str]:
    """The Threadfin warnings of the last :func:`call_path` / :func:`run` / :func:`render` call (R issues them as R
    warnings)."""
    return list(_LAST_WARNINGS)


def _function(path: Any) -> Any:
    import threadfin

    obj: Any = threadfin
    for name in [path] if isinstance(path, str) else list(path):
        obj = getattr(obj, str(name))
    return obj


def call_path(path: Any, args: Any = None, kwargs: Any = None) -> Any:
    """Call the Threadfin function at ``path`` (e.g. ``["tl", "clone_profiles"]``), R values coerced."""
    fn = _function(path)
    args, kwargs = coerce_call(fn, tuple(args or ()), dict(kwargs or {}))
    return call_recording(fn, args, kwargs)


def attr(obj: Any, name: str, default: Any = None) -> Any:
    """``getattr(obj, name, default)`` (R reaches attributes without reticulate converting them)."""
    return getattr(obj, name, default)


def method(obj: Any, name: str, args: Any = None, kwargs: Any = None) -> Any:
    """Call ``obj.<name>(*args, **kwargs)``."""
    return getattr(obj, name)(*tuple(args or ()), **dict(kwargs or {}))


def item(obj: Any, key: Any) -> Any:
    """``obj[key]``."""
    return obj[key]


def keys(obj: Any) -> list[str]:
    """The keys of a mapping, as strings."""
    return [str(k) for k in dict(obj or {})]


def write_text(path: str, text: str) -> str:
    """Write ``text`` to ``path`` as UTF-8, as ``pathlib.Path.write_text`` does."""
    from pathlib import Path

    Path(path).write_text(str(text), encoding="utf-8")
    return path


# --------------------------------------------------------------------------------------------- objects from R data


def _index(values: Any) -> pd.Index:
    return pd.Index(np.atleast_1d(np.asarray(values, dtype=object)).astype(str), dtype=object)


def _names(values: Any) -> list[str]:
    return [] if values is None else [str(v) for v in np.atleast_1d(np.asarray(values, dtype=object))]


def frame_from_r(frame: Any, *, index: Any = None, logical: Any = None, set_index: Any = None) -> pd.DataFrame:
    """A pandas table from an R data frame converted by reticulate, with R's missing values restored.

    R integer and logical columns with ``NA`` arrive as doubles (the R side converts them, because reticulate turns
    those ``NA`` into -2147483648 / ``False``); ``logical`` names the columns to turn back into pandas' nullable
    booleans. ``index`` sets the row labels (strings); ``set_index`` names columns that become the index (the
    index of a table that came from pandas).
    """
    df = frame if isinstance(frame, pd.DataFrame) else pd.DataFrame(frame)
    df = df.copy()
    for c in _names(logical):
        v = pd.to_numeric(df[c], errors="coerce").to_numpy(dtype=float)
        df[c] = pd.array([None if np.isnan(x) else bool(x) for x in v], dtype="boolean")
    if index is not None:
        df.index = _index(index)
    cols = [c for c in _names(set_index) if c in df.columns]
    if cols:
        df = df.set_index(cols)
        if len(cols) == 1 and df.index.name == "index":
            df.index.name = None
    return df


def contigs_from_r(columns: Any, names: Any) -> pd.DataFrame:
    """A receptor contig table from R as pandas reads the same table from a file.

    R gives every column as text, written as it would stand in a file (``None`` for a missing value). The table is
    written to CSV in memory and read back with ``pandas.read_csv(..., low_memory=False)``, the call with which
    Threadfin reads ``filtered_contig_annotations.csv`` and AIRR files, so that column types, missing values and
    booleans come out exactly as when Threadfin reads the file itself.
    """
    names = _names(names)
    cols = {n: np.asarray(list(c) if isinstance(c, list | tuple) else np.atleast_1d(c), dtype=object)
            for n, c in zip(names, columns)}
    text = pd.DataFrame(cols, columns=names).to_csv(index=False)
    return pd.read_csv(io.StringIO(text), low_memory=False)


def anndata_from_r(
    x: Any,
    i: Any,
    p: Any,
    n_obs: int,
    n_vars: int,
    obs_names: Any,
    var_names: Any,
    obs: pd.DataFrame | None = None,
    obsm: dict[str, Any] | None = None,
) -> Any:
    """AnnData of raw counts from the slots of an R ``dgCMatrix`` (genes x cells), obs metadata and embeddings.

    The compressed columns of the genes x cells matrix are the compressed rows of its cells x genes transpose, so
    ``(x, i, p)`` become a CSR matrix of shape ``(n_obs, n_vars)`` without reordering a value. Whole-number counts
    below 2**24 are stored as float32 (exact; the dtype of scanpy's readers), others as float64. ``obsm`` maps keys
    to cells x dimensions matrices (Seurat reductions, reducedDims).
    """
    import anndata as ad
    from scipy import sparse

    data = np.asarray(x, dtype=np.float64)
    if data.size and np.all(np.mod(data, 1) == 0) and np.abs(data).max() < 2**24:
        data = data.astype(np.float32)
    elif not data.size:
        data = data.astype(np.float32)
    m = sparse.csr_matrix(
        (data, np.asarray(i, dtype=np.int32), np.asarray(p, dtype=np.int32)), shape=(int(n_obs), int(n_vars))
    )
    names = _index(obs_names)
    if obs is None or (isinstance(obs, pd.DataFrame) and obs.shape[1] == 0):
        obs = pd.DataFrame(index=names)
    else:
        obs = obs.copy()
        obs.index = names
    var = pd.DataFrame(index=_index(var_names))
    a = ad.AnnData(X=m, obs=obs, var=var)
    for key, value in dict(obsm or {}).items():
        a.obsm[str(key)] = np.asarray(value, dtype=float).reshape(int(n_obs), -1)
    return a


def read_h5ad(path: str) -> Any:
    """``anndata.read_h5ad(path)``."""
    import anndata as ad

    return ad.read_h5ad(path)


# ----------------------------------------------------------------------------------------------- the one call


def result_from_adata(adata: Any) -> Any:
    """``threadfin.ThreadfinResult.from_adata(adata)``, with Threadfin's warnings kept for R."""
    import threadfin

    _LAST_WARNINGS.clear()
    return call_recording(threadfin.ThreadfinResult.from_adata, (adata,), {})


def run_columns(result: Any) -> pd.DataFrame:
    """The cell-level results of :func:`threadfin.run`, indexed by cell: the ``obs`` columns it writes (its
    docstring): the clone ids, the ``bcr_*`` receptor columns when ``bcr`` was given, ``clone_programme`` and
    ``clone_programme_assigned``, in the order of ``obs``."""
    adata = result.adata
    clone_key = str(result.settings.get("clone_key") or "clone_id")
    cols = [c for c in adata.obs.columns
            if str(c) in (clone_key, "clone_programme", "clone_programme_assigned") or str(c).startswith("bcr_")]
    out = adata.obs.loc[:, cols].copy()
    out.index = pd.Index(adata.obs_names.astype(str), name="cell", dtype=object)
    return out


def series_from_r(values: Any, names: Any, name: Any = None, attrs: Any = None) -> pd.Series:
    """A pandas Series (object values, string index) from an R named vector, with ``attrs`` restored (a clone-id
    vector from :func:`threadfin.define_clones_dandelion` carries its settings there)."""
    vals = [None if (v is None or (isinstance(v, float) and np.isnan(v))) else v
            for v in (list(values) if isinstance(values, list | tuple) else np.atleast_1d(values).tolist())]
    s = pd.Series(vals, index=_index(names), name=None if name is None else str(name), dtype=object)
    s.index.name = "cell_id"
    s.attrs.update(dict(attrs or {}))
    return s


# ----------------------------------------------------------------------------------------------- tables for R


def _default_index(index: pd.Index) -> bool:
    return isinstance(index, pd.RangeIndex) and index.start == 0 and index.step == 1 and index.name is None


def _text(v: Any) -> Any:
    if v is None or (isinstance(v, float) and np.isnan(v)) or v is pd.NA or v is pd.NaT:
        return None
    if isinstance(v, str):
        return v
    if isinstance(v, list | tuple | set | frozenset | np.ndarray | pd.Series):
        return ";".join(str(x) for x in v)
    return str(v)


def frame_for_r(df: Any, index: bool | None = None) -> tuple[pd.DataFrame | None, list[str], list[str]]:
    """A copy of ``df`` in dtypes reticulate converts faithfully, the columns R should read as logical, and the
    columns that held the index.

    Strings (any pandas string dtype or object column) become object columns of ``str`` with ``None`` for missing
    values (lists are joined with ``";"``); nullable integers become int64 (float64 with NaN when values are
    missing); nullable booleans become bool (float64 1 / 0 / NaN when values are missing, listed as logical);
    categoricals keep their codes with string categories. The index becomes the first column(s) when it carries
    labels (``index=None``: anything but a default ``RangeIndex``).
    """
    if df is None:
        return None, [], []
    if isinstance(df, pd.Series):
        df = df.to_frame(name=df.name if df.name is not None else "value")
    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df)
    keep = (not _default_index(df.index)) if index is None else bool(index)
    index_cols: list[str] = []
    if keep and isinstance(df.index, pd.MultiIndex):
        index_cols = [str(n) if n is not None else f"level_{k}" for k, n in enumerate(df.index.names)]
        df = df.reset_index()
    elif keep:
        name = str(df.index.name) if df.index.name is not None else "index"
        if name in df.columns:
            name = f"{name}_index"
        df = df.reset_index(names=name)
        index_cols = [name]
    out: dict[str, Any] = {}
    logical: list[str] = []
    for c in df.columns:
        s = df[c]
        name = str(c)
        dt = s.dtype
        if isinstance(dt, pd.CategoricalDtype):
            cats = [str(x) for x in dt.categories]
            out[name] = pd.Categorical.from_codes(s.cat.codes.to_numpy(), categories=pd.Index(cats, dtype=object),
                                                  ordered=bool(dt.ordered))
        elif pd.api.types.is_bool_dtype(dt) and not isinstance(dt, pd.BooleanDtype):
            out[name] = s.to_numpy(dtype=bool)
        elif isinstance(dt, pd.BooleanDtype):
            if s.isna().any():
                out[name] = s.astype("float64").to_numpy(dtype=float, na_value=np.nan)
                logical.append(name)
            else:
                out[name] = s.to_numpy(dtype=bool)
        elif pd.api.types.is_integer_dtype(dt):
            if s.isna().any():
                out[name] = s.to_numpy(dtype=float, na_value=np.nan)
            else:
                out[name] = _r_integers(s.to_numpy(dtype=np.int64))
        elif pd.api.types.is_float_dtype(dt):
            out[name] = s.to_numpy(dtype=float, na_value=np.nan)
        elif pd.api.types.is_datetime64_any_dtype(dt) or pd.api.types.is_timedelta64_dtype(dt):
            out[name] = s.astype(str).where(s.notna(), None).to_numpy(dtype=object)
        elif pd.api.types.is_complex_dtype(dt):
            out[name] = np.asarray([_text(v) for v in s.to_numpy(dtype=object)], dtype=object)
        else:  # strings, objects, mixed
            vals = s.to_numpy(dtype=object, na_value=None)
            if all(isinstance(v, bool | np.bool_) for v in vals if v is not None) and any(v is not None for v in vals):
                if any(v is None for v in vals):
                    out[name] = np.asarray([np.nan if v is None else float(v) for v in vals], dtype=float)
                    logical.append(name)
                else:
                    out[name] = vals.astype(bool)
            else:
                out[name] = np.asarray([_text(v) for v in vals], dtype=object)
    return pd.DataFrame(out, index=pd.RangeIndex(len(df))), logical, index_cols


def _r_integers(a: np.ndarray) -> np.ndarray:
    """Whole numbers as int32 when they fit, which R reads as integers (int64 arrives as double)."""
    if a.size == 0 or (a.min() >= -(2**31 - 1) and a.max() <= 2**31 - 1):
        return a.astype(np.int32)
    return a


def _scalar(v: Any) -> Any:
    if isinstance(v, np.generic):
        return v.item()
    if v is pd.NA or v is pd.NaT:
        return None
    return v


def _plain(v: Any) -> bool:
    """A value made only of strings, numbers, booleans and None (also in lists, tuples and dicts)."""
    if v is None or isinstance(v, str | bool | int | float | np.generic):
        return True
    if isinstance(v, list | tuple):
        return all(_plain(x) for x in v)
    if isinstance(v, dict):
        return all(_plain(x) for x in v.values())
    return False


def _attrs_for_r(attrs: dict) -> dict:
    """The entries of ``DataFrame.attrs`` / ``Series.attrs`` that are plain values: settings, counts, flags. Arrays,
    tables and other objects kept there (the clone bootstrap of a gene table) stay in Python, where the functions
    that need them find them."""
    return {str(k): value_for_r(v) for k, v in dict(attrs or {}).items() if _plain(v)}


def value_for_r(obj: Any) -> Any:
    """Any Threadfin result as plain values R converts faithfully.

    dicts stay dicts, tables become ``{"__frame__": table, "__logical__": [...], "__index__": [...],
    "__attrs__": {...}}`` (see :func:`frame_for_r`), a Series becomes ``{"__series__": values, "__names__": index,
    "__attrs__": {...}}``, arrays and lists or tuples of numbers, strings or booleans become NumPy arrays (R
    vectors), lists of records (dicts with the same keys) become tables, NumPy scalars become Python scalars and
    anything else is left as it is (a Python object in R).
    """
    if obj is None:
        return None
    if isinstance(obj, pd.DataFrame):
        frame, logical, index_cols = frame_for_r(obj)
        return {"__frame__": frame, "__logical__": logical, "__index__": index_cols, "__attrs__": _attrs_for_r(obj.attrs)}
    if isinstance(obj, pd.Series):
        vals = obj.to_numpy(dtype=object, na_value=None)
        if pd.api.types.is_numeric_dtype(obj.dtype) and not pd.api.types.is_bool_dtype(obj.dtype):
            vals = obj.to_numpy(dtype=float, na_value=np.nan)
        elif pd.api.types.is_bool_dtype(obj.dtype) and not obj.isna().any():
            vals = obj.to_numpy(dtype=bool)
        else:
            vals = np.asarray([_text(v) for v in vals], dtype=object)
        return {"__series__": vals, "__names__": np.asarray([str(i) for i in obj.index], dtype=object),
                "__name__": None if obj.name is None else str(obj.name), "__attrs__": _attrs_for_r(obj.attrs)}
    if isinstance(obj, dict):
        return {str(k): value_for_r(v) for k, v in obj.items()}
    if isinstance(obj, np.ndarray):
        if obj.ndim == 0:                   # a single value, e.g. read back from an h5ad file
            return value_for_r(obj.item())
        if obj.dtype.kind in "iu":
            return _r_integers(obj.astype(np.int64)) if obj.dtype.itemsize >= 4 else obj
        if obj.dtype.kind in "bf":
            return obj
        return np.asarray([_text(v) for v in obj.ravel()], dtype=object).reshape(obj.shape)
    if isinstance(obj, list | tuple):
        items = [_scalar(v) for v in obj]
        if items and all(isinstance(v, dict) for v in items) and len({tuple(v) for v in items}) == 1:
            return value_for_r(pd.DataFrame(items))
        if all(isinstance(v, bool) for v in items) and items:
            return np.asarray(items, dtype=bool)
        if all(isinstance(v, int | float) and not isinstance(v, bool) or v is None for v in items) and items:
            if any(v is None for v in items):
                return np.asarray([np.nan if v is None else float(v) for v in items], dtype=float)
            if all(isinstance(v, int) for v in items):
                return _r_integers(np.asarray(items, dtype=np.int64))
            return np.asarray(items, dtype=float)
        if all(isinstance(v, str) or v is None for v in items) and items:
            return np.asarray(items, dtype=object)
        return [value_for_r(v) for v in items]
    return _scalar(obj)


def frame_value(obj: Any) -> Any:
    """A table as :func:`value_for_r` gives it; a Series becomes a one-column table first."""
    if obj is None:
        return None
    if isinstance(obj, pd.Series):
        obj = obj.to_frame(name=obj.name if obj.name is not None else "value")
    if not isinstance(obj, pd.DataFrame):
        obj = pd.DataFrame(obj)
    return value_for_r(obj)


def matrix_for_r(m: Any) -> dict[str, Any]:
    """Slots of the genes x cells ``dgCMatrix`` of a cells x genes matrix (CSR arrays = CSC of the transpose)."""
    from scipy import sparse

    csr = sparse.csr_matrix(m)
    csr.sum_duplicates()
    csr.sort_indices()
    return {
        "x": csr.data.astype(np.float64),
        "i": csr.indices.astype(np.int32),
        "p": csr.indptr.astype(np.int32),
        "n_obs": int(csr.shape[0]),
        "n_vars": int(csr.shape[1]),
    }


# ------------------------------------------------------------------------------------------------------- figures


def _figure(out: Any) -> Any:
    import matplotlib.pyplot as plt
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure, SubFigure

    if isinstance(out, Figure):
        return out, True
    if isinstance(out, SubFigure):
        return out.figure, True
    if isinstance(out, Axes):
        return out.figure, False
    if isinstance(out, np.ndarray | list | tuple) and len(out):
        first = np.asarray(out, dtype=object).ravel()[0]
        if isinstance(first, Axes | Figure):
            return _figure(first)
    if out is None and plt.get_fignums():
        return plt.gcf(), True
    raise TypeError(
        f"the plotting function returned {type(out).__name__}, not a matplotlib Figure or Axes; call it directly "
        "(tf$pl$...) when it writes files itself"
    )


def render(fn: Any, args: Any, kwargs: Any, save: Any = None, png: Any = None, dpi: Any = None,
           target: Any = None) -> dict[str, Any]:
    """Call a Threadfin plotting function and close every figure it opened.

    ``save`` is the file to write. Without ``dpi`` it is passed to the function itself (every ``threadfin.pl``
    function and ``ThreadfinResult.plot`` take it), so that the file is the one Threadfin writes, in the format of
    its extension; with ``dpi`` the figure is saved here at that resolution, with the same matplotlib settings.
    ``png``, when given, is a PNG of the same figure for the R graphics device, at ``dpi`` or, with ``target`` =
    (width, height) in pixels, at the resolution at which the figure fills the device.
    """
    if "matplotlib.pyplot" not in sys.modules:
        import matplotlib

        matplotlib.use("Agg")
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt

    before = set(plt.get_fignums())
    fn = getattr(fn, "__wrapped__", fn)
    kwargs = dict(kwargs or {})
    if save is not None and dpi is None:
        kwargs["save"] = str(save)
    args, kwargs = coerce_call(fn, tuple(args or ()), kwargs)
    try:
        fig, _ = _figure(call_recording(fn, args, kwargs))
        background = mcolors.to_hex(fig.get_facecolor(), keep_alpha=False)
        width, height = (float(v) for v in fig.get_size_inches())
        if save is not None and dpi is not None:
            fig.savefig(str(save), dpi=float(dpi))
        if png is not None:
            res = dpi
            if res is None and target is not None:
                tw, th = (float(v) for v in target)
                res = max(20.0, min(tw / width, th / height))
            fig.savefig(str(png), format="png", dpi="figure" if res is None else float(res),
                        facecolor=fig.get_facecolor(), bbox_inches="tight")
    finally:
        for n in set(plt.get_fignums()) - before:
            plt.close(n)
    return {"save": None if save is None else str(save), "png": None if png is None else str(png),
            "background": background, "width_in": width, "height_in": height}


def read_png(path: str) -> np.ndarray:
    """RGBA pixels of a PNG as uint8 (height x width x 4), for R without the png package."""
    import matplotlib.image as mpimg

    img = mpimg.imread(path)
    if img.dtype != np.uint8:
        img = np.clip(np.round(img * 255), 0, 255).astype(np.uint8)
    if img.ndim == 2:
        img = np.stack([img] * 3 + [np.full_like(img, 255)], axis=-1)
    elif img.shape[2] == 3:
        img = np.concatenate([img, np.full(img.shape[:2] + (1,), 255, dtype=np.uint8)], axis=-1)
    return np.ascontiguousarray(img)
