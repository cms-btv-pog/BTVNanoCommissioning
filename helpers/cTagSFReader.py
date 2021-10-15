import uproot
import numpy as np
import coffea
from coffea import hist, processor
def getSF(flav, CvL, CvB, file="DeepCSV_ctagSF_MiniAOD94X_2017_pTincl.root",syst=""):    

    #print "Reading from %s. Systematic: %s." %(file,syst)
    
    f = uproot.open(file)
    
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
