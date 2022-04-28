from BTVNanoCommissioning.helpers.definitions import definitions
import hist as Hist
def histogrammer(workflow):
    _hist_dict={}
    ## Common variables
    flav_axis = Hist.axis.IntCategory([0, 1, 4, 5, 6], name="flav", label="Genflavour")
    syst_axis = Hist.axis.StrCategory([], name="syst", growth=True)
    pt_axis = Hist.axis.Regular(50, 0, 300, name="pt", label=" $p_{T}$ [GeV]")
    softlpt_axis = Hist.axis.Regular(25, 0, 25, name="pt", label=" $p_{T}$ [GeV]")
    mass_axis = Hist.axis.Regular(50, 0, 300, name="mass", label=" $p_{T}$ [GeV]")
    eta_axis = Hist.axis.Regular(25, -2.5, 2.5, name="eta", label=" $\eta$")
    phi_axis = Hist.axis.Regular(30, -3, 3, name="phi", label="$\phi$")
    mt_axis = Hist.axis.Regular(30, 0, 300, name="mt", label=" $m_{T}$ [GeV]")
    iso_axis = Hist.axis.Regular(30, 0, 0.15, name="pfRelIso03_all", label="Rel. Iso")
    softliso_axis = Hist.axis.Regular(20, 0.2, 4.2, name="pfRelIso03_all", label="Rel. Iso")
    dr_axis = Hist.axis.Regular(20, 0, 4, name="dr", label="$\Delta$R")
    sliso_axis = Hist.axis.Regular(30, 0.2, 4., name="pfRelIso03_all", label="Rel. Iso")
    dxy_axis = Hist.axis.Regular(40, -0.05, 0.05, name="dxy", label="d_{xy}")
    dz_axis = Hist.axis.Regular(40, 0, 0.1, name="dz", label="d_{z}")
    sip3d_axis = Hist.axis.Regular(20, 0, 0.2, name="sip3d", label="SIP 3D")
    ptratio_axis = Hist.axis.Regular(50, 0, 1, name="ratio", label="ratio")
    n_axis = Hist.axis.IntCategory([0, 1, 2, 3, 4, 5], name="n", label="N obj")
    
    
    ### Workflow specific 
    if "ttdilep_sf" == workflow :
        obj_list = ["mu","ele"]
        for i in range(2):
            obj_list.append(f"jet{i}")
            _hist_dict[f"dr_mujet{i}"] =  Hist.Hist(flav_axis,dr_axis,Hist.storage.Weight())
        for i in ["mu","ele"]:
            _hist_dict[f"{i}_pfRelIso04_all"] = Hist.Hist(iso_axis,Hist.storage.Weight())
            _hist_dict[f"{i}_dxy"] = Hist.Hist(dxy_axis,Hist.storage.Weight())
            _hist_dict[f"{i}_dz"] = Hist.Hist(dz_axis,Hist.storage.Weight())
    if "ttsemilep_sf" == workflow:
        obj_list = ["mu","MET"]
        for i in range(4):
            obj_list.append(f"jet{i}")
            _hist_dict[f"dr_mujet{i}"] =  Hist.Hist(flav_axis,dr_axis,Hist.storage.Weight())
        for i in ["mu"]:
            _hist_dict[f"{i}_pfRelIso04_all"] = Hist.Hist(iso_axis,Hist.storage.Weight())
            _hist_dict[f"{i}_dxy"] = Hist.Hist(dxy_axis,Hist.storage.Weight())
            _hist_dict[f"{i}_dz"] = Hist.Hist(dz_axis,Hist.storage.Weight())
    elif "ctag_ttdilep_sf" in workflow:
        obj_list = ["hl","sl","soft_l","MET","z","lmujet"]
        _hist_dict["z_mass"] = Hist.Hist(Hist.axis.Regular(50, 50, 100, name="mass", label=" $p_{T}$ [GeV]"),Hist.storage.Weight())
        
        _hist_dict["dr_lmujetsmu"] =  Hist.Hist(flav_axis,dr_axis,Hist.storage.Weight())
        _hist_dict["dr_lmujethmu"] =  Hist.Hist(flav_axis,dr_axis,Hist.storage.Weight())
        _hist_dict["dr_lmusmu"] =  Hist.Hist(dr_axis,Hist.storage.Weight())
        for i in ["hl","sl","soft_l"]:
            if i == "soft_l":_hist_dict[f"soft_l_pfRelIso04_all"] = Hist.Hist(softliso_axis,Hist.storage.Weight())
            else :_hist_dict[f"{i}_pfRelIso04_all"] = Hist.Hist(iso_axis,Hist.storage.Weight())
            _hist_dict[f"{i}_ptratio"] = Hist.Hist(flav_axis,ptratio_axis,Hist.storage.Weight())
            _hist_dict[f"{i}_dxy"] = Hist.Hist(dxy_axis,Hist.storage.Weight())
            _hist_dict[f"{i}_dz"] = Hist.Hist(dz_axis,Hist.storage.Weight())
    elif "ctag_ttsemilep_sf" in workflow or "Wc_sf" in workflow:
        obj_list = ["hl","soft_l","MET","z","w","mujet"]
        _hist_dict["z_mass"] = Hist.Hist(Hist.axis.Regular(50, 50, 100, name="mass", label=" $p_{T}$ [GeV]"),Hist.storage.Weight())
        _hist_dict["w_mass"] = Hist.Hist(Hist.axis.Regular(50, 50, 100, name="mass", label=" $p_{T}$ [GeV]"),Hist.storage.Weight())
        
        _hist_dict["dr_lmujetsmu"] =  Hist.Hist(flav_axis,dr_axis,Hist.storage.Weight())
        _hist_dict["dr_lmujethmu"] =  Hist.Hist(flav_axis,dr_axis,Hist.storage.Weight())
        _hist_dict["dr_lmusmu"] =  Hist.Hist(dr_axis,Hist.storage.Weight())
        for i in ["hl","soft_l"]:
            if i == "soft_l":_hist_dict[f"soft_l_pfRelIso04_all"] = Hist.Hist(softliso_axis,Hist.storage.Weight())
            else :_hist_dict[f"{i}_pfRelIso04_all"] = Hist.Hist(iso_axis,Hist.storage.Weight())
            _hist_dict[f"{i}_ptratio"] = Hist.Hist(flav_axis,ptratio_axis,Hist.storage.Weight())
            _hist_dict[f"{i}_dxy"] = Hist.Hist(dxy_axis,Hist.storage.Weight())
            _hist_dict[f"{i}_dz"] = Hist.Hist(dz_axis,Hist.storage.Weight())
    elif "DY_sf" in workflow:
        obj_list = ["posl","negl","z","jet"]
        _hist_dict["z_mass"] = Hist.Hist(Hist.axis.Regular(50, 50, 100, name="mass", label=" $p_{T}$ [GeV]"),Hist.storage.Weight())
        _hist_dict["dr_mumu"] =  Hist.Hist(dr_axis,Hist.storage.Weight())
        for i in ["posl","negl"]:
            _hist_dict[f"{i}_pfRelIso04_all"] = Hist.Hist(iso_axis,Hist.storage.Weight())
            _hist_dict[f"{i}_dxy"] = Hist.Hist(dxy_axis,Hist.storage.Weight())
            _hist_dict[f"{i}_dz"] = Hist.Hist(dz_axis,Hist.storage.Weight())
    ### Common kinematic variables

    _hist_dict["njet"] = Hist.Hist(n_axis,Hist.storage.Weight())
    for obj in obj_list:
        if "jet" in obj : 
            _hist_dict[f"{obj}_pt"] = Hist.Hist(flav_axis,pt_axis,Hist.storage.Weight())
            _hist_dict[f"{obj}_eta"] = Hist.Hist(flav_axis,eta_axis,Hist.storage.Weight())
            _hist_dict[f"{obj}_phi"] = Hist.Hist(flav_axis,phi_axis,Hist.storage.Weight())
            _hist_dict[f"{obj}_mass"] = Hist.Hist(flav_axis,mass_axis,Hist.storage.Weight())
        elif obj == "MET":
            _hist_dict[f"{obj}_pt"] = Hist.Hist(pt_axis,Hist.storage.Weight())
            _hist_dict[f"{obj}_phi"] = Hist.Hist(phi_axis,Hist.storage.Weight())
        else:
            if obj == "soft_l":_hist_dict["soft_l_pt"] = Hist.Hist(softlpt_axis,Hist.storage.Weight())
            else:_hist_dict[f"{obj}_pt"] = Hist.Hist(pt_axis,Hist.storage.Weight())
            _hist_dict[f"{obj}_phi"] = Hist.Hist(phi_axis,Hist.storage.Weight())
            _hist_dict[f"{obj}_eta"] = Hist.Hist(eta_axis,Hist.storage.Weight())
            
    ### Btag variables
    bininfo = definitions()
    for d in bininfo.keys():
        ranges = bininfo[d]["manual_ranges"]
        binning = bininfo[d]["bins"]
        if ranges[1] is None:
            ranges[1] = 0.0
        if ranges[0] is None:
            ranges[0] = -0.5
        _hist_dict[d] = Hist.Hist(
        flav_axis,
        Hist.axis.Regular(binning, ranges[0], ranges[1], name=d, label=d),
        Hist.storage.Weight())
    ### discriminators        
    disc_list = ["btagDeepB","btagDeepC","btagDeepFlavB","btagDeepFlavC","btagDeepCvL","btagDeepCvB","btagDeepFlavCvL","btagDeepFlavCvB"]
    for disc in disc_list:
        nlep=1
        if "ttdilep_sf" in workflow: nlep = 2
        elif "ttsemilep_sf" in workflow: nlep = 4
        for i in range(nlep):
            _hist_dict[f"{disc}_{i}"] = Hist.Hist(
            flav_axis,
            syst_axis,
            Hist.axis.Regular(30,  -0.2, 1, name="discr", label="discr"),
            Hist.storage.Weight())
    return  _hist_dict
        
     
