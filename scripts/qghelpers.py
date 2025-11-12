import numpy as np
import hist as Hist
import itertools
import copy
from coffea.util import load

flavs = ["ud", "s", "c", "b", "g", "other"]

def load_inputs(coffea_files):
    input_hists = {}
    available_vars = set()
    objects_in_hists = set()
    for f in coffea_files:
        if ".coffea" not in f:
            print(f"Skipping non-coffea file {f}")
            continue

        fdict = load(f)
        for sample in fdict:
            if sample not in input_hists:
                input_hists[sample] = fdict[sample]
                for var in fdict[sample]:
                    if not var.endswith("_pteta"):
                        continue
                    available_vars.add(var)
                    objects_in_hists.add(var.split("_Var")[0].replace("Obj",""))
            else:
                for var in fdict[sample]:
                    if not var.endswith("_pteta"):
                        continue
                    available_vars.add(var)
                    objects_in_hists.add(var.split("_Var")[0].replace("Obj",""))
                    if var not in input_hists[sample]:
                        input_hists[sample][var] = fdict[sample][var]
                    else:
                        input_hists[sample][var] += fdict[sample][var]

    # Get pT and eta distributions as projections
    for sample in input_hists:
        for obj in objects_in_hists:
            pt_var = f"Obj{obj}_Varpt"
            eta_var = f"Obj{obj}_Vareta"
            for var in input_hists[sample]:
                if f"Obj{obj}_Var" not in var or not var.endswith("_pteta"):
                    continue
                if "flav" in input_hists[sample][var].axes.name:
                    pr = ["syst", "flav", "pt", "eta"]
                else:
                    pr = ["syst", "pt", "eta"]
                input_hists[sample][pt_var] = input_hists[sample][var].project(*pr)
                available_vars.add(pt_var)
                input_hists[sample][eta_var] = input_hists[sample][var].project(*pr)
                available_vars.add(eta_var)
                break

    return input_hists, available_vars

def get_slices(range_str, bincoords=False):
    slices = []
    parts = range_str.split(",")
    for r in range(len(parts)-1):
        start, end = parts[r], parts[r+1]
        if bincoords:
            slices.append(slice(int(start), int(end)))
        else:
            slices.append(slice(complex(0, float(start)), complex(0, float(end))))
    return slices

def stack_hists(hdict, vars_to_plot, regions= [{"pt": slice(None), "eta": slice(None)}]):
    stacked = {}
    var_regions = list(itertools.product(vars_to_plot, regions))

    for var, region in var_regions:
        region_str = "_".join([f"{k.upper()}{v.start}-{v.stop}" for k, v in region.items()])
        print(f"Stacking variable {var} in region {region_str}")
        stacked[f"{var}_{region_str}"] = {"data": {}, "mc_flav": {}, "mc_sample": {}}
        for sample in hdict:
            if var not in hdict[sample]:
                print(f"Variable {var} not found in sample {sample}, skipping")
                continue

            h = hdict[sample][var][region]
            if "_pteta" in var:
                summed_h = h[{"pt": slice(0, len(h.axes["pt"]), Hist.rebin(len(h.axes["pt"]))),
                    "eta": slice(0, len(h.axes["eta"]), Hist.rebin(len(h.axes["eta"])))
                    }]
                summed_h = summed_h[{"pt": 0, "eta": 0}]
            elif "_Varpt" in var:
                summed_h = h[{"eta": slice(0, len(h.axes["eta"]), Hist.rebin(len(h.axes["eta"])))}]
                summed_h = summed_h[{"eta": 0}]
            elif "_Vareta" in var:
                summed_h = h[{"pt": slice(0, len(h.axes["pt"]), Hist.rebin(len(h.axes["pt"])))}]
                summed_h = summed_h[{"pt": 0}]

            # Sum over flavor axis for data
            if "Run" in sample:
                h_dt = summed_h[{"flav": sum}] if "flav" in summed_h.axes.name else summed_h
                stacked[f"{var}_{region_str}"]["data"][sample] = h_dt
                continue

            h_mc_sample = summed_h[{"flav": sum}] if "flav" in summed_h.axes.name else summed_h
            if sample not in stacked[f"{var}_{region_str}"]["mc_flav"]:
                stacked[f"{var}_{region_str}"]["mc_sample"][sample] = h_mc_sample
            else:
                stacked[f"{var}_{region_str}"]["mc_sample"][sample] += h_mc_sample

            for i, flav in enumerate(flavs):
                if "flav" in summed_h.axes.name:
                    h_mc_flav = summed_h[{"flav": i}]
                else:
                    h_mc_flav = summed_h
                if flav not in stacked[f"{var}_{region_str}"]["mc_flav"]:
                    stacked[f"{var}_{region_str}"]["mc_flav"][flav] = h_mc_flav
                else:
                    stacked[f"{var}_{region_str}"]["mc_flav"][flav] += h_mc_flav

    return stacked

def plot_roc(hists, systname="nominal", var="Var"):
    mc_total = sum([hists["mc_flav"][flav] for flav in hists["mc_flav"]])
    # Gluon jets as gluon and other flav
    ghist = hists["mc_flav"]["g"]# + hists["mc_flav"]["other"]
    ghist = ghist[{"syst": systname}].values()
    # Quark jets as udscb
    flavs = ["ud", "s",]# "c", "b"]
    qhist = sum([hists["mc_flav"][flav] for flav in flavs])
    qhist = qhist[{"syst": systname}].values()

    # Normalize to 1
    norm_factor = (qhist + ghist)
    tpr = []
    fpr = []
    for i, (g, q) in enumerate(zip(ghist, qhist)):
        if i == 0:
            tpr += [1.0]
            fpr += [1.0]
            continue
        if i == len(ghist) - 1:
            tpr += [0.0]
            fpr += [0.0]
            continue
        tpr_temp = np.sum(qhist[i:]) / np.sum(qhist)
        fpr_temp = np.sum(ghist[i:]) / np.sum(ghist)
        if fpr_temp == fpr[-1]:
            continue
        tpr += [tpr_temp]
        fpr += [fpr_temp]

    from sklearn.metrics import auc
    # Sort by fpr
    tpr, fpr = zip(*sorted(zip(tpr, fpr)))
    auc = np.trapz(tpr, fpr)

    import matplotlib.pyplot as plt
    plt.plot(fpr, tpr, label="ROC curve")
    plt.plot([0, 1], [0, 1], linestyle="--", color="gray", label="Random guess")
    plt.text(0.6, 0.2, f"AUC = {auc:.3f}", fontsize=12, bbox=dict(facecolor='white', alpha=0.5))
    plt.xlabel("False Positive Rate (Gluon efficiency)")
    plt.ylabel("True Positive Rate (Quark efficiency)")
    plt.title("Quark vs Gluon Jet Discrimination ROC")
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.grid()
    plt.legend()
    plt.savefig(f"qg_roc_curve_{var}.png")
    plt.close()
    

    
