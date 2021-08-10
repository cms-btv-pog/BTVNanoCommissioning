from .DummyTagger1 import DummyTagger1
from .DummyTagger2 import DummyTagger2

taggers = {}

taggers["DummyTagger1"] = DummyTagger1
taggers["DummyTagger2"] = DummyTagger2

__all__ = ["taggers"]
