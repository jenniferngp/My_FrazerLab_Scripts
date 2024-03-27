def plot_gc(data_file):
    '''
    Replot the Picard output as png file to put into the html
    '''
    # Load data
    data = pd.read_table(data_file, comment="#")

    # Plot the data
    fig = plt.figure()
    ax = fig.add_subplot(111)

    plt.xlim((0, 100))

    lin1 = ax.plot(data['GC'], data['NORMALIZED_COVERAGE'],
                   label='Normalized coverage', color='r')
    ax.set_ylabel('Normalized coverage')

    ax2 = ax.twinx()
    lin2 = ax2.plot(data['GC'], data['MEAN_BASE_QUALITY'],
                    label='Mean base quality at GC%', color='b')
    ax2.set_ylabel('Mean base quality at GC%')

    ax3 = ax.twinx()
    lin3 = ax3.plot(data['GC'], data['WINDOWS']/np.sum(data['WINDOWS']),
                    label='Windows at GC%', color='g')
    ax3.get_yaxis().set_visible(False)

    lns = lin1 + lin2 + lin3
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc='best')

    # plot_img = BytesIO()
    # fig.savefig(plot_img, format='png')
    prefix = data_file.rstrip('.gc.txt')
    plot_png = prefix + '.gc_plot.png'
    fig.savefig(plot_png, format='png')
    
    print("{} saved!".format(plot_png))
    

import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
data_file = sys.argv[1]

print("Data file:" + data_file)

plot_gc(data_file)