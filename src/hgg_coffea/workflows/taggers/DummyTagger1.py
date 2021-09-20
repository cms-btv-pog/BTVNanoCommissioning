import awkward
import numpy


class DummyTagger1:
    def __init__(self) -> None:
        pass

    @property
    def name(self) -> str:
        return "DummyTagger1"

    # lower priority is better
    # first decimal point is category within tag (in case for untagged)
    @property
    def priority(self) -> int:
        return 20

    def __call__(self, events: awkward.Array) -> awkward.Array:
        ##-- Baseline example for subcategorization:
        """
        ##-- Derive new column from logic 
        # tag per diphoton, not per event
#         events.diphotons
        counts = ak.num(events.diphoton.pt, axis=1) ##-- Number of entries per row. (N diphotons per row)  
        cats = numpy.random.randint(low=0, high=4, size = len(awkward.flatten(events.diphoton.pt)) ##-- Number of diphotons 
        cats = ak.unflatten(cats, counts) ##-- Back to size of events ##-- Want to select highest pT diphoton object ##-- Zero out diphotons not being used. 
        return (self.priority + cats) * ....
        return (self.priority + 1) * awkward.ones_like(events.diphotons.pt, dtype=numpy.int32) ##-- returning a column of all 20s. This is categorizing the events. Giving entire tagger outcome. 
        """ 
        return (
            self.priority * awkward.ones_like(events.diphotons.pt, dtype=numpy.int32),
            {},
        )