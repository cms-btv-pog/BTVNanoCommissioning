"""
20 September 2021
Abraham Tishelman-Charny

The purpose of this tagger is to tag HH->WWgg events.
"""

import awkward as ak
import numpy


class HHWWggTagger:
    def __init__(self) -> None:
        pass

    @property
    def name(self) -> str:
        return "HHWWggTagger"

    # lower priority is better
    # first decimal point is category within tag (in case for untagged)
    @property
    def priority(self) -> int:
        return 40

    # Lepton selections
    def electron_selection(self, electrons: ak.Array) -> ak.Array:

        return electrons[(electrons.pt > self.min_pt_electron)]

    def muon_selection(self, muons: ak.Array) -> ak.Array:

        return muons[(muons.pt > self.min_pt_muon)]

    # Jet selections
    def jet_selection(self, jets: ak.Array) -> ak.Array:

        return jets[
            (jets.pt > self.min_pt_jet) & (abs(jets.eta) < self.max_abs_eta_jet)
        ]

    def GetCategory(self, ievt: int) -> int:

        evt_nLeptons = self.nLeptons[ievt]
        evt_nJets = self.nJets[ievt]

        if evt_nLeptons == 1:  # Semileptonic
            cat = 0
        elif evt_nLeptons == 0 and evt_nJets >= 4:  # Fullyhadronic
            cat = 1
        elif evt_nLeptons >= 2:  # Fullyleptonic
            cat = 2
        else:  # Untagged
            cat = 3

        return cat

    def __call__(self, events: ak.Array) -> ak.Array:

        """
        Pass the lepton and jet selections for a given tagger through the metaconditions?
        """

        # Electron selections
        self.min_pt_electron = 10.0

        # Muon selections
        self.min_pt_muon = 10.0

        # Jet selections
        self.min_pt_jet = 30.0
        self.max_abs_eta_jet = 2.4

        # Get base lepton collections for selections
        electrons = events.Electron
        muons = events.Muon

        electrons = self.electron_selection(electrons)
        muons = self.muon_selection(muons)

        # Get jets for jet selection
        jets = events.Jet
        jets = self.jet_selection(jets)

        nElectrons = ak.num(electrons, axis=1)
        nMuons = ak.num(muons, axis=1)
        nLeptons = numpy.add(nElectrons, nMuons)
        nJets = ak.num(jets, axis=1)

        self.nLeptons = nLeptons
        self.nJets = nJets

        nDiphotons = ak.num(
            events.diphotons.pt, axis=1
        )  # Number of entries per row. (N diphotons per row)
        ievts_by_dipho = ak.flatten(
            ak.Array([nDipho * [evt_i] for evt_i, nDipho in enumerate(nDiphotons)])
        )
        cat_vals = ak.Array(map(self.GetCategory, ievts_by_dipho))
        cats = ak.unflatten(cat_vals, nDiphotons)  # Back to size of events.
        cats_by_diphoEvt = self.priority + cats

        return (cats_by_diphoEvt, {})
