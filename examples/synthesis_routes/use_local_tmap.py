import sys
import os

# Add your module paths
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT_DIR = os.path.abspath(os.path.join(_SCRIPT_DIR, '../..'))

# Add the src directory to the Python path
sys.path.insert(0, os.path.join(_ROOT_DIR, 'src'))

# Now directly import the local modules
from tmap.core.synthesis_tree import SynthesisTree, SynthesisTreeNode
from tmap.layout_generators.ted_layout_generator import TEDLayoutGenerator
from tmap.helpers.ted_utils import compute_ted

# Make these available in the global namespace
# (No need to import tmap module since we're using direct imports)