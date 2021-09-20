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
    
    def GetCategory(self, ievt) -> int:
        
        evt_nLeptons = self.nLeptons[ievt] 
        evt_nJets = self.nJets[ievt]
        
        if(evt_nLeptons == 0):
            cat = 0 
        elif(evt_nLeptons == 1 and evt_nJets >= 4):
            cat = 1 
        elif(evt_nLeptons >= 2):
            cat = 2
        else:
            cat = 3 
            
        return cat 

    def __call__(self, events: ak.Array) -> ak.Array:
        
        electrons = events["electrons"]
        muons = events["muons"]
        jets = events["jets"]
        
        nElectrons = ak.num(electrons, axis=1)
        nMuons = ak.num(muons, axis=1)
        nLeptons = numpy.add(nElectrons, nMuons)
        nJets = ak.num(jets, axis=1)
        
        self.nLeptons = nLeptons
        self.nJets = nJets
        
        ievts = numpy.array([i for i in range(len(events))])

        nDiphotons = ak.num(events.diphotons.pt, axis=1) ##-- Number of entries per row. (N diphotons per row) 
        ievts_by_dipho = ak.flatten(ak.Array([nDipho*[evt_i] for evt_i, nDipho in enumerate(nDiphotons)]))
        cat_vals = ak.Array(map(self.GetCategory, ievts_by_dipho))
        print("cat_vals:",cat_vals)

#         cats = numpy.random.randint(low=0, high=4, size = len(ak.flatten(events.diphotons.pt))) ##-- Randomly assign categories. The length of this array is the number of diphotons 
        cats = ak.unflatten(cat_vals, nDiphotons) ##-- Back to size of events. 
        
        cats_by_diphoEvt = self.priority + cats 
        print("cats_by_diphoEvt:",cats_by_diphoEvt)
        
        return (cats_by_diphoEvt, {})
