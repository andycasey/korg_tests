from grok.transitions import Transitions
from grok.photospheres import Photosphere
from grok.photospheres.interpolator import PhotosphereInterpolator

DO_NOT_SCALE_ABUNDANCES = ("H", "He") # This is the same as what SME does. Disentangling it might be a little tricky, so we do the same.