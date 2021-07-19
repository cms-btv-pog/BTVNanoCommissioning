import awkward as ak
import numpy as np


class DummyTagger1:
    def __init__(self):
        pass

    def name(self):
        return self.__class__

    # lower priority is better
    def priority(self):
        return 20

    def __call__(self, events):
        x = np.random.uniform(low=-1.0, high=1.0, size=len(events))
        return x > 0
