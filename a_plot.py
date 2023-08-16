from matplotlib.pyplot import *

def a_subplot(nrows=1,ncols=2,dpi=150,figsize=(10, 6),pad=3.0):
    figure
    fig, axs = subplots(nrows,ncols,dpi=dpi,figsize=figsize)
    fig.tight_layout(pad=pad)
    return axs
