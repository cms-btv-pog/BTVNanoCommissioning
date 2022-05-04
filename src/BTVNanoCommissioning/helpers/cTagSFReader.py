import uproot
import numpy as np
import awkward as ak
import importlib
import coffea
from coffea import hist, processor

# FIXME - use coffea dense_lookup instead
def getSF(flav, CvL, CvB, file="DeepCSV_ctagSF_MiniAOD94X_2017_pTincl.root",syst=""):    
    if not isinstance(flav, np.ndarray):
        flav = ak.to_numpy(flav)
    if not isinstance(CvL, np.ndarray):
        CvL = ak.to_numpy(CvL)
    if not isinstance(CvB, np.ndarray):
        CvB = ak.to_numpy(CvB)

    _btag_path = "BTVNanoCommissioning.data.BTV.Rereco17_94X"
    with importlib.resources.path(_btag_path, file) as filename:
        f = uproot.open(filename)
    
    if syst == "" or syst == "central": systsuff = ""
    else: systsuff = '_'+syst
    
    bins = f['SFl_hist'+systsuff].to_numpy()[-1]
    SFd = {
        0: f['SFl_hist'+systsuff].to_numpy(flow=True)[0],
        4: f['SFc_hist'+systsuff].to_numpy(flow=True)[0],
        5: f['SFb_hist'+systsuff].to_numpy(flow=True)[0],
    }

    icvl = np.digitize(np.clip(CvL, -1, 0.999999), [-1]+list(bins)) - 1
    icvb = np.digitize(np.clip(CvB, -1, 0.999999), [-1]+list(bins)) - 1

    SFarr = np.empty(len(flav))
    for key in set(flav):
        ix = np.where(np.array(flav) == key)
        SFarr[ix] = SFd[key][icvl[ix], icvb[ix]]
        
    if "Stat" in syst: return np.absolute(SFarr - getSF(flav,CvL,CvB,file))
    return SFarr
