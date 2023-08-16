from numpy import *
from matplotlib.pyplot import *

import os
cwd = os.getcwd()

from scipy import constants as cst
from scipy import integrate
from scipy import optimize
from scipy import signal
from scipy import special
from scipy import interpolate
from matplotlib.cm import ScalarMappable

fon = 10;
matplotlib.rc('xtick', labelsize=fon) 
matplotlib.rc('ytick', labelsize=fon)
rcParams['axes.titlesize'] = fon
matplotlib.rc('legend',fontsize=fon-4)
from IPython.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))
