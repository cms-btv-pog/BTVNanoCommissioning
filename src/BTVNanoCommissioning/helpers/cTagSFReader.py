import uproot
import numpy as np
import importlib.resources
from coffea.lookup_tools.dense_lookup import dense_lookup
import awkward as ak


def getSF(flav, CvL, CvB, file="DeepCSV_ctagSF_MiniAOD94X_2017_pTincl.root", syst=""):
    # _btag_path = "BTVNanoCommissioning.data.BTV.Rereco17_94X"
    with importlib.resources.path(
        file[: file.rfind("/")].replace("/", "."), file[file.rfind("/") + 1 :]
    ) as filename:
        f = uproot.open(filename)

    if syst == "" or syst == "central":
        systsuff = ""
    else:
        systsuff = "_" + syst

    SFd = np.array(
        [
            f["SFl_hist" + systsuff].to_numpy()[0],
            f["SFc_hist" + systsuff].to_numpy()[0],
            f["SFb_hist" + systsuff].to_numpy()[0],
        ]
    )
    bins = [
        f["SFb_hist" + systsuff].to_numpy()[-1],
        f["SFb_hist" + systsuff].to_numpy()[-2],
        np.array([0, 1, 2, 3]),
    ]
    efflookup = dense_lookup(SFd, bins)
    flav = np.where(flav == 4, 1, flav)
    flav = np.where(flav == 5, 2, flav)
    CvL = ak.to_numpy(CvL)
    CvB = ak.to_numpy(CvB)
    flav = ak.to_numpy(flav)
    SFarr = efflookup(CvL, CvB, flav)
    if "Stat" in syst:
        return np.absolute(SFarr - getSF(flav, CvL, CvB, file))
    return SFarr
