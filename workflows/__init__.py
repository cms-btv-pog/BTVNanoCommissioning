from .dystudies import DYStudiesProcessor
from . import taggers

workflows = {}

workflows['dystudies'] = DYStudiesProcessor

__all__ = ['workflows', 'taggers']
