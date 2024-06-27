#############################################################
# REQUIRED PACAKGES
#############################################################
from ipylar.version import *
import os

#############################################################
# INITIALIZATION
#############################################################
print("Loading iPyLAR version",version)

#############################################################
# DATA RETRIEVAL
#############################################################
#Root directory
try:
    FILE=__file__
    ROOTDIR=os.path.abspath(os.path.dirname(FILE))
except:
    FILE=""
    ROOTDIR=os.path.abspath('')

def get_data_path(filename):
    return os.path.join(ROOTDIR, 'data', filename)

def list_data():
    return os.listdir(os.path.join(ROOTDIR, 'data'))

#############################################################
# GLOBAL CONSTANTS
#############################################################
class Const(object):
    # Time-zone
    UTCL=5*3600
    
    # Time constants
    hours=3600
    days=86400
    years=365.25*days

    # Physical constants
    
#############################################################
# LOAD MODULES
#############################################################
from ipylar.basins import *
