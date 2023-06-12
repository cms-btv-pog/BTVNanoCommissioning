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
    return out


def num(ar):
    return ak.num(ak.fill_none(ar[~ak.is_none(ar)], 0), axis=0)


## Based on https://github.com/CoffeaTeam/coffea/discussions/735
def _is_rootcompat(a):
    """Is it a flat or 1-d jagged array?"""
    t = ak.type(a)
    if isinstance(t, ak._ext.ArrayType):
        if isinstance(t.type, ak._ext.PrimitiveType):
            return True
        if isinstance(t.type, ak._ext.ListType) and isinstance(
            t.type.type, ak._ext.PrimitiveType
        ):
            return True
    return False


def uproot_writeable(events, include=["events", "run", "luminosityBlock"]):
    ev = {}
    include = np.array(include)
    no_filter = False

    if len(include) == 1 and include[0] == "*":
        no_filter = False
    for bname in events.fields:
        if not events[bname].fields:
            if not no_filter and bname not in include:
                continue
            ev[bname] = ak.fill_none(
                ak.packed(ak.without_parameters(events[bname])), -99
            )
        else:
            b_nest = {}
            no_filter_nest = False
            if all(np.char.startswith(include, bname) == False):
                continue
            include_nest = [
                i[i.find(bname) + len(bname) + 1 :]
                for i in include
                if i.startswith(bname)
            ]

            if len(include_nest) == 1 and include_nest[0] == "*":
                no_filter_nest = True

            if not no_filter_nest:
                mask_wildcard = np.char.find(include_nest, "*") != -1
                include_nest = np.char.replace(include_nest, "*", "")
            for n in events[bname].fields:
                ## make selections to the filter case, keep cross-ref ("Idx")
                if (
                    not no_filter_nest
                    and all(np.char.find(n, include_nest) == -1)
                    and "Idx" not in n
                    and "Flavor" not in n
                ):
                    continue
                if not _is_rootcompat(events[bname][n]) and ak.num(
                    events[bname][n], axis=0
                ) != len(flatten(events[bname][n])):
                    continue

                b_nest[n] = ak.fill_none(
                    ak.packed(ak.without_parameters(events[bname][n])), -99
                )
            ev[bname] = ak.zip(b_nest)
    return ev
