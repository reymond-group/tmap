from .mstmap import *

import os

_ROOT = os.path.abspath(os.path.dirname(__file__))
def get_asset(path):
    return os.path.join(_ROOT, 'assets', path)