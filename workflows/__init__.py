from .dystudies import DYStudiesProcessor
from .taggers import taggers

workflows = {}

workflows["dystudies"] = DYStudiesProcessor

__all__ = ["workflows", "taggers"]
