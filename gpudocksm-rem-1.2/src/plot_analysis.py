import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

def plotTotalEne(values, pdf_path):
    """line plot of the total energy
    """
    ofn = pdf_path
    ener = values

    fig, ax = plt.subplots()
    ax.plot(ener)
    ax.set_xlabel('step', fontsize=14)
    ax.set_ylabel('total energy', fontsize=14)

    print "line plot of total energy\t", ofn
    plt.savefig(ofn)
    
def histogram(values, pdf_path, normed=1, num_bins=100):
    """plot the histogram of the data
    
    Arguments:
    - `values`: data
    - `pdf_path`: path to save the out put
    - `normed`: flag, 1 means to normalize the bin height, 0 does not
    """
    mu = values.mean()
    sigma = values.std()
    x = values
    
    fig, ax = plt.subplots()

    # the histogram of the data
    num_bins = num_bins
    n, bins, patches = ax.hist(x, num_bins, normed=normed, facecolor='green', alpha=0.5)

    # add a 'best fit' line
    y = mlab.normpdf(bins, mu, sigma)
    ax.plot(bins, y, 'r--')
    ax.set_xlabel('Totol energy')
    ax.set_ylabel('Probability')
    ax.set_title(r'$\mu$=' + str(mu) + ', ' + r'$\sigma$=' + str(sigma))

    # Tweak spacing to prevent clipping of ylabel
    plt.subplots_adjust(left=0.15)
    print "histgram plot of total energy\t", pdf_path
    fig.savefig(pdf_path)








