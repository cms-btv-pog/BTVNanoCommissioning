from typing import List, Optional, Tuple

import awkward
import numpy
import xgboost


def calculate_diphoton_mva(
    mva: Tuple[Optional[xgboost.Booster], List[str]],
    diphotons: awkward.Array,
    events: awkward.Array,
) -> awkward.Array:
    if mva[0] is None:
        return diphotons
    diphoton_mva = mva[0]

    var_order = mva[1]

    bdt_vars = {}

    bdt_vars["dipho_leadIDMVA"] = diphotons.pho_lead.mvaID
    bdt_vars["dipho_subleadIDMVA"] = diphotons.pho_sublead.mvaID
    bdt_vars["dipho_leadEta"] = diphotons.pho_lead.eta
    bdt_vars["dipho_subleadEta"] = diphotons.pho_sublead.eta
    bdt_vars["dipho_lead_ptoM"] = diphotons.pho_lead.pt / diphotons.mass
    bdt_vars["dipho_sublead_ptoM"] = diphotons.pho_sublead.pt / diphotons.mass

    def calc_displacement(
        photons: awkward.Array, events: awkward.Array
    ) -> awkward.Array:
        x = photons.x_calo - events.PV.x
        y = photons.y_calo - events.PV.y
        z = photons.z_calo - events.PV.z
        return awkward.zip({"x": x, "y": y, "z": z}, with_name="Vector3D")

    v_lead = calc_displacement(diphotons.pho_lead, events)
    v_sublead = calc_displacement(diphotons.pho_sublead, events)

    p_lead = v_lead.unit() * diphotons.pho_lead.energyRaw
    p_lead["energy"] = diphotons.pho_lead.energyRaw
    p_lead = awkward.with_name(p_lead, "Momentum4D")
    p_sublead = v_sublead.unit() * diphotons.pho_sublead.energyRaw
    p_sublead["energy"] = diphotons.pho_sublead.energyRaw
    p_sublead = awkward.with_name(p_sublead, "Momentum4D")

    sech_lead = 1.0 / numpy.cosh(p_lead.eta)
    sech_sublead = 1.0 / numpy.cosh(p_sublead.eta)
    tanh_lead = numpy.cos(p_lead.theta)
    tanh_sublead = numpy.cos(p_sublead.theta)

    cos_dphi = numpy.cos(p_lead.deltaphi(p_sublead))

    numerator_lead = sech_lead * (
        sech_lead * tanh_sublead - tanh_lead * sech_sublead * cos_dphi
    )
    numerator_sublead = sech_sublead * (
        sech_sublead * tanh_lead - tanh_sublead * sech_lead * cos_dphi
    )

    denominator = 1.0 - tanh_lead * tanh_sublead - sech_lead * sech_sublead * cos_dphi

    add_reso = (
        0.5
        * (-numpy.sqrt(2.0) * events.BeamSpot.sigmaZ / denominator)
        * (numerator_lead / p_lead.mag + numerator_sublead / p_sublead.mag)
    )

    dEnorm_lead = diphotons.pho_lead.energyErr / diphotons.pho_lead.energy
    dEnorm_sublead = diphotons.pho_sublead.energyErr / diphotons.pho_sublead.energy

    sigma_m = 0.5 * numpy.sqrt(dEnorm_lead ** 2 + dEnorm_sublead ** 2)
    sigma_wv = numpy.sqrt(add_reso ** 2 + sigma_m ** 2)

    vtx_prob = awkward.full_like(sigma_m, 0.999)  # !!!! placeholder !!!!

    bdt_vars["CosPhi"] = cos_dphi
    bdt_vars["vtxprob"] = vtx_prob
    bdt_vars["sigmarv"] = sigma_m
    bdt_vars["sigmawv"] = sigma_wv

    counts = awkward.num(diphotons, axis=-1)
    bdt_inputs = numpy.column_stack(
        [awkward.to_numpy(awkward.flatten(bdt_vars[name])) for name in var_order]
    )
    tempmatrix = xgboost.DMatrix(bdt_inputs, feature_names=var_order)
    scores = diphoton_mva.predict(tempmatrix)

    for var, arr in bdt_vars.items():
        if "dipho" not in var:
            diphotons[var] = arr

    diphotons["bdt_score"] = awkward.unflatten(scores, counts)

    return diphotons
