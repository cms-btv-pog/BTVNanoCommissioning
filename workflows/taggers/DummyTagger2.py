import awkward as ak
import numpy as np


class DummyTagger2:
    def __init__(self):
        pass

    def name(self):
        return self.__class__

    # lower priority is better
    def priority(self):
        return 10

    def __call__(self, events):
        x = np.random.exponential(size=len(events))
        return x > 2
