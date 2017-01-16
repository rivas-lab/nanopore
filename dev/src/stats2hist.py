#!/usr/bin/env python2

_README_ = '''
-------------------------------------------------------------------------
plot read length distribution
-------------------------------------------------------------------------
'''

import pandas as pd
import matplotlib, argparse
matplotlib.use('agg')
from matplotlib import pyplot as plt

def make_hist(x, nbins = 20, title = None, xlabel = None, ylabel = None, filename = None):
    '''
    This function generates histogram of a vector x and save to file
    Inputs:
      x: data vector
      nbins:    number of bins
      title:    title of the plot
      xlabel:   label on x-axis
      ylabel:   label on y-axis
      filename: name of the image file (if given, save to file)
    Returns:
      matlab plot object
    Side effect:
      save an image file if filename is given
    '''
    
    import matplotlib
    matplotlib.use('agg')
    from matplotlib import pyplot as plt
    
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(1, 1, 1)
    ax.hist(x, nbins)
    
    if(xlabel != None):
        ax.set_xlabel(xlabel)
    if(ylabel != None):
        ax.set_ylabel(ylabel)
    if(title != None):
        ax.set_title(title)
    if(filename != None):
        fig.savefig(filename)
    
    return(fig)

def stats2hist_main(data_f, nbins):
    df = pd.read_csv(data_f, sep = ' ', names=['chr', 'MAPQ', 'len'])
    make_hist(x = df.ix[:, 2],
              nbins = nbins,
              title='Read length distribution of PacBio NA12878',
              xlabel='read length [nt]', ylabel = 'count', 
              filename = data_f + '.hist.png')
    return(0)


def main():    
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=_README_)
    
    parser.add_argument('-b', metavar='B', type=int,
                        default = 100,
                        help='number of bins (default: 100)')
    
    parser.add_argument('i', metavar='I',                     
                        help='input stats file name')        
    args = parser.parse_args()
    
    return(stats2hist_main(data_f = args.i, nbins = args.b))   


if __name__ == "__main__":
    main()
