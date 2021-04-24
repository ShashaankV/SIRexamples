

from matplotlib import rc
from matplotlib import rcParams
import matplotlib.pyplot as plt
plt.ion()
import seaborn as sns

# rcParams['axes.facecolor'] = 'white'

def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

def plot_style(fs=14):
    font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : fs}
    rc('font', **font)
    rcParams['axes.spines.right'] = False
    rcParams['axes.spines.left'] = False
    rcParams['axes.spines.bottom'] = False
    rcParams['axes.spines.top'] = False

def plot_w_ptlab(x,y,ptlab,xlab=None,ylab=None,fs=10,ax=None):
    ax.plot(x,y,'.')
    for i in range(len(ptlab)):
        ax.annotate(ptlab[i],xy=(x[i],y[i]),fontsize=fs)
    if xlab:
        ax.set_xlabel(xlab)
    if ylab:
        ax.set_ylabel(ylab)
    