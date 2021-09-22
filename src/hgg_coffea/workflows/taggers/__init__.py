from typing import Any, Dict

from .DummyTagger1 import DummyTagger1
from .DummyTagger2 import DummyTagger2
from .HHWWggTagger import HHWWggTagger

taggers: Dict[str, Any] = {}

taggers["DummyTagger1"] = DummyTagger1
taggers["DummyTagger2"] = DummyTagger2
taggers["HHWWggTagger"] = HHWWggTagger

__all__ = ["taggers"]
