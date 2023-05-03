import numpy as np
import matplotlib.pyplot as plt
import os
import gc
from tqdm import tqdm

def plot_cell_fractions(input_dir, output_dir, runs=20):
    """
    Function that plots fractions of cells with each ecDNA type
    
    Arguments:
    input_dir (str) : path to directory containing output of C++ simulations
    output_dir (str): path to directory where output plots will be saved
    runs (int)      : number of simulation runs (default 20)
    
    Returns:
    The function returns None if run successfully. Saves plots to output_dir
    """
    # only select files from input_dir that contain cellFractions in the filename
    files = [f for f in os.listdir(input_dir) if (os.path.isfile(os.path.join(input_dir, f)) and \
                                               f.split('_')[0]=='cellFractions')]
    # sort files according to reverse date (earliest first)
    files.sort(key=lambda f: os.path.getmtime(os.path.join(input_dir,f)))
    # create output_dir if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    for f in tqdm(files):
        # plot data
        data_fracs = np.loadtxt(os.path.join(input_dir, f),
                                dtype=float)
        data_fracs = np.reshape(data_fracs, [runs,-1,5])
        
        
        for d in data_fracs:
            plt.plot(d[:,1], d[:,2], '--', c='tab:blue', linewidth=0.2)
            plt.plot(d[:,1], d[:,3], '--', c='tab:orange', linewidth=0.2)
            plt.plot(d[:,1], d[:,4], '--', c='tab:green', linewidth=0.2)
        means = np.mean(data_fracs, axis=0)
        plt.plot(means[:,1], means[:,2], '-', c='tab:blue', linewidth=2, label='a')
        plt.plot(means[:,1], means[:,3], '-', c='tab:orange', linewidth=2, label='b')
        plt.plot(means[:,1], means[:,4], '-', c='tab:green', linewidth=2, label='neutral')
        plt.legend(loc=1)
        plt.xscale('log')
        
        plt.savefig(os.path.join(
            output_dir, f.split('.')[0])+'.png')
        plt.close()

        # explicitly clear memory
        del data_fracs
        gc.collect()
    
    return None

