import awkward as ak
import numpy as np


def flatten(ar):  # flatten awkward into a 1d array to hist
    return ak.flatten(ar, axis=None)


def normalize(val, cut):
    if cut is None:
        ar = ak.to_numpy(ak.fill_none(val, np.nan))
        return ar
    else:
        ar = ak.to_numpy(ak.fill_none(val[cut], np.nan))
        return ar


def update(events, collections):
    """Return a shallow copy of events array with some collections swapped out"""
    out = events
    for name, value in collections.items():
        out = ak.with_field(out, value, name)
        if (
            "float" not in str(ak.type(value))
            and "int" not in str(ak.type(value))
            and "bool" not in str(ak.type(value))
        ):
            if name == "Jet":
                out.Jet["pt"] = value.pt
            elif name == "MET":
                out.MET["pt"] = value.pt
                out.MET["phi"] = value.phi
    return out


def num(ar):
    return ak.num(ak.fill_none(ar[~ak.is_none(ar)], 0), axis=0)
