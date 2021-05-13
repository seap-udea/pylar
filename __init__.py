import numpy as np
from scipy.integrate import odeint,solve_ivp
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import rc
rc("animation",html="jshtml")
from matplotlib.patches import Rectangle
import warnings
warnings.filterwarnings("ignore")