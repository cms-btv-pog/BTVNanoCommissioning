from hgg_coffea.workflows.dystudies import DYStudiesProcessor
from hgg_coffea.workflows.taggers import taggers

workflows = {}

workflows["dystudies"] = DYStudiesProcessor

__all__ = ["workflows", "taggers", "DYStudiesProcessor"]
