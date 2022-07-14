import sys
import os
import math
import shared
import mcell as m

MODEL_PATH = os.path.dirname(os.path.abspath(__file__))

n_glu = 100 #number of glutamate per vesicles


# ---- model parameters ----

# declare all items from parameter_overrides as variables
for parameter_name, value in shared.parameter_overrides.items():
    setattr(sys.modules[__name__], parameter_name, value)

# auxiliary function used to determine whether a parameter was defined
def not_defined(parameter_name):
    return parameter_name not in globals()

# ---- simulation setup ----

if not_defined('ITERATIONS'):
    ITERATIONS = 1e4
