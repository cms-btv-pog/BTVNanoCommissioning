import gzip
import warnings
from typing import Optional

import xgboost


def _get_gzip(fname: str) -> bytearray:
    return bytearray(gzip.open(fname, "rb").read())


_magics = {b"\x1f\x8b": _get_gzip}


def load_bdt(fname: str) -> Optional[xgboost.Booster]:
    try:
        bdt = xgboost.Booster()
        with open(fname, "rb") as f:
            magic = f.read(2)
            opener = _magics.get(magic, lambda x: x)
        bdt.load_model(opener(fname))
    except xgboost.core.XGBoostError as xgberr:
        warnings.warn(repr(xgberr))
        bdt = None
    except FileNotFoundError as fnferr:
        warnings.warn(repr(fnferr))
        bdt = None
    return bdt
