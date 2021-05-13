import numpy as np
from scipy.integrate import odeint,solve_ivp
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import rc
rc("animation",html="jshtml")
from matplotlib.patches import Rectangle
import warnings
warnings.filterwarnings("ignore")

#############################################
#GRAPHICAL PARAMETERS
#############################################
#Properties
#Size of the box
width=10
#Vertical offset text
voff=0.2
#Margin boxes 
fmar=1
#Separation box
fsep=2
#Height box
hrep=3
#Margin vertical channels
mver=0.2

#Decorative properties
fs=12
#Fluid in repositories
wc='c'
#Fluid in channels
tc='c'
#Color lines of channels
cc='k'
#Width of walls
lw=2

#Size repositories and channels
uL=hrep
uA=hrep
unq=1
uQ=1
uP=1
uE=1