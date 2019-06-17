import os
import sys
import pickle
import math

import umap
import numpy as np
import pandas as pd
import scipy.stats as ss

from timeit import default_timer as timer
from mnist import MNIST
from progress.bar import Bar
from progress.spinner import Spinner
from mhfp.encoder import MHFPEncoder
from tmap import tmap
from operator import itemgetter
from faerun import Faerun

import matplotlib.pyplot as plt
from matplotlib import colors


labels, values = pickle.load(open('pdb.pickle', 'rb'))
values.append(np.array(values[2]) / np.array(values[0]))

plt.rcParams['mathtext.fontset'] = 'dejavuserif'
VERY_SMALL_SIZE = 6
SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIG_SIZE = 12
plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=SMALL_SIZE)
plt.rc('axes', labelsize=SMALL_SIZE)
plt.rc('xtick', labelsize=VERY_SMALL_SIZE)
plt.rc('ytick', labelsize=VERY_SMALL_SIZE)
plt.rc('legend', fontsize=SMALL_SIZE)
plt.rc('figure', titlesize=BIG_SIZE)

_, axarr = plt.subplots(ncols=2, figsize=(8.27 / 1.0, 8.27 / 11.0), dpi=600)

N1, bins1, patches1 = axarr[0].hist(values[0], bins=5000)
N2, bins2, patches2 = axarr[1].hist(values[7], bins=5000)

total = 0
for thisN, thispatch, thisbin in zip(N1, patches1, bins1):
    color = plt.cm.rainbow(total / len(values[0]))
    thispatch.set_facecolor(color)
    total += thisN

total = 0
for thisN, thispatch, thisbin in zip(N2, patches2, bins2):
    color = plt.cm.rainbow(total / len(values[7]))
    thispatch.set_facecolor(color)
    total += thisN

axarr[0].get_yaxis().set_visible(False)
axarr[1].get_yaxis().set_visible(False)

axarr[0].set_xlabel('Size (Heavy Atom Count)')
axarr[1].set_xlabel('Fraction of Negative Charges')

axarr[0].set_yscale('log', nonposy='clip')
axarr[1].set_yscale('log', nonposy='clip')

axarr[0].set_xscale('log', nonposx='clip')
axarr[1].set_xscale('log', nonposx='clip')

# axarr[0].set_xlim(0, 25000)

plt.tight_layout()
plt.savefig('labels_pdb.png')