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

    def __call__(self, events: ak.Array) -> ak.Array:
        
        counts = ak.num(events.diphotons.pt, axis=1) ##-- Number of entries per row. (N diphotons per row) 
        cats = numpy.random.randint(low=0, high=4, size = len(ak.flatten(events.diphotons.pt))) ##-- Randomly assign categories. The length of this array is the number of diphotons 
        cats = ak.unflatten(cats, counts) ##-- Back to size of events. 
        
        default_cats = self.priority * ak.ones_like(events.diphotons.pt, dtype=numpy.int32)
        updated_cats = self.priority + cats 
        
        return (updated_cats, {})

        """
        ##-- Derive new column from logic 
        # tag per diphoton, not per event
#         events.diphotons
        counts = ak.num(events.diphoton.pt, axis=1) ##-- Number of entries per row. (N diphotons per row)  
        cats = numpy.random.randint(low=0, high=4, size = len(ak.flatten(events.diphoton.pt)) ##-- Number of diphotons 
        cats = ak.unflatten(cats, counts) ##-- Back to size of events ##-- Want to select highest pT diphoton object ##-- Zero out diphotons not being used. 
        return (self.priority + cats) * ....
        return (self.priority + 1) * ak.ones_like(events.diphotons.pt, dtype=numpy.int32) ##-- returning a column of all 20s. This is categorizing the events. Giving entire tagger outcome. 
        """ 
        return (
            self.priority * ak.ones_like(events.diphotons.pt, dtype=numpy.int32),
            {},
        )