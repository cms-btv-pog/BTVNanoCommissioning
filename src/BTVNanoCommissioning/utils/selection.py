import awkward as ak


## Jet pu ID not exist in Winter22Run3 sample
def jet_id(events, campaign):
    if campaign == "Rereco17_94X":
        jetmask = (
            (events.Jet.pt > 20)
            & (abs(events.Jet.eta) <= 2.5)
            & (events.Jet.jetId >= 5)
            & ((events.Jet.pt > 50) | (events.Jet.puId >= 7))
        )
    else:
        jetmask = (
            (events.Jet.pt > 20)
            & (abs(events.Jet.eta) <= 2.5)
            & (events.Jet.jetId >= 5)
        )
    return jetmask


## FIXME: Electron cutbased Id & MVA ID not exist in Winter22Run3 sample
def ele_cuttightid(events, campaign):
    elemask = (
        (abs(events.Electron.eta) < 1.4442)
        | ((abs(events.Electron.eta) < 2.5) & (abs(events.Electron.eta) > 1.566))
    ) & (events.Electron.cutBased > 3)
    return elemask


def ele_mvatightid(events, campaign):
    elemask = (
        (abs(events.Electron.eta) < 1.4442)
        | ((abs(events.Electron.eta) < 2.5) & (abs(events.Electron.eta) > 1.566))
    ) & (events.Electron.mvaIso_WP80 > 0.5)

    return elemask


def softmu_mask(events, campaign):
    softmumask = (
        (events.Muon.pt < 25)
        & (abs(events.Muon.eta) < 2.4)
        & (events.Muon.tightId > 0.5)
        & (events.Muon.pfRelIso04_all > 0.2)
        & (events.Muon.jetIdx != -1)
    )

    return softmumask


def mu_idiso(events, campaign):
    mumask = (
        (abs(events.Muon.eta) < 2.4)
        & (events.Muon.tightId > 0.5)
        & (events.Muon.pfRelIso04_all <= 0.15)
    )
    return mumask


def btag_mu_idiso(events, campaign):
    mumask = (
        (abs(events.Muon.eta) < 2.4)
        & (events.Muon.tightId > 0.5)
        & (events.Muon.pfRelIso04_all < 0.12)
    )
    return mumask


def jet_cut(events, campaign):
    multijetmask = (
        (abs(events.Jet.eta) < 2.4) & (events.Jet.pt > 180) & (events.Jet.jetId >= 5)
    )
    return multijetmask


btag_wp_dict = {
    "Summer22Run3": {
        "DeepFlav": {
            "b": {
                "No": 0.0,
                "L": 0.0583,
                "M": 0.3086,
                "T": 0.7183,
                "XT": 0.8111,
                "XXT": 0.9512,
            },
            "c": {
                "No": [0.0, 0.0],
                "L": [0.042, 0.208],  # CvL, then CvB
                "M": [0.108, 0.299],
                "T": [0.303, 0.243],
            },
        },
        "RobustParTAK4": {
            "b": {
                "No": 0.0,
                "L": 0.0849,
                "M": 0.4319,
                "T": 0.8482,
                "XT": 0.9151,
                "XXT": 0.9874,
            },
            "c": {
                "No": [0.0, 0.0],
                "L": [0.039, 0.068],
                "M": [0.117, 0.130],
                "T": [0.360, 0.095],
            },
        },
        "PNet": {
            "b": {
                "No": 0.0,
                "L": 0.047,
                "M": 0.245,
                "T": 0.6734,
                "XT": 0.7862,
                "XXT": 0.961,
            },
            "c": {
                "No": [0.0, 0.0],
                "L": [0.054, 0.181],
                "M": [0.160, 0.306],
                "T": [0.492, 0.259],
            },
        },
    },
    "Summer22EERun3": {
        "DeepFlav": {
            "b": {
                "No": 0.0,
                "L": 0.0614,
                "M": 0.3196,
                "T": 0.73,
                "XT": 0.8184,
                "XXT": 0.9542,
            },
            "c": {
                "No": [0.0, 0.0],
                "L": [0.042, 0.206],
                "M": [0.108, 0.298],
                "T": [0.305, 0.241],
            },
        },
        "RobustParTAK4": {
            "b": {
                "No": 0.0,
                "L": 0.0897,
                "M": 0.451,
                "T": 0.8604,
                "XT": 0.9234,
                "XXT": 0.9893,
            },
            "c": {
                "No": [0.0, 0.0],
                "L": [0.039, 0.067],
                "M": [0.117, 0.128],
                "T": [0.358, 0.095],
            },
        },
        "PNet": {
            "b": {
                "No": 0.0,
                "L": 0.0499,
                "M": 0.2605,
                "T": 0.6915,
                "XT": 0.8033,
                "XXT": 0.9664,
            },
            "c": {
                "No": [0.0, 0.0],
                "L": [0.054, 0.182],
                "M": [0.160, 0.304],
                "T": [0.491, 0.258],
            },
        },
    },
}


def btag_wp(jets, campaign, tagger, borc, wp):
    WP = btag_wp_dict[campaign]
    if borc == "b":
        jet_mask = jets[f"btag{tagger}B"] > WP[tagger]["b"][wp]
    else:
        jet_mask = (jets[f"btag{tagger}CvB"] > WP[tagger]["c"][wp][1]) & (
            jets[f"btag{tagger}CvL"] > WP[tagger]["c"][wp][0]
        )
    return jet_mask
