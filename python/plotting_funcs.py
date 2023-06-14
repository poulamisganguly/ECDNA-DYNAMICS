import numpy as np
import matplotlib.pyplot as plt
import os
import gc
from tqdm import tqdm

def plot_cell_fractions(input_dir, output_dir, runs, num_cells, only_means):
    """
    Function that plots fractions of cells with each ecDNA type
    
    Arguments:
    input_dir (str)  : path to directory containing output of C++ simulations
    output_dir (str) : path to directory where output plots will be saved
    runs (int)       : number of simulation runs (default 20)
    num_cells (int)  : number of cells in final population (default 100)
    only_means (bool): plot only mean trajectory (default True)
    
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
        data_fracs = np.reshape(data_fracs, [runs,num_cells,-1])
        fig1, ax1 = plt.subplots()
        fig2, ax2 = plt.subplots() 
        if not only_means:
            for d in data_fracs:
                for i in range(2,data_fracs.shape[2]-1):
                    ax1.plot(d[:,1], d[:,i], '--', linewidth=0.2)
                ax1.plot(d[:,1], d[:,-1], '--', linewidth=0.2, c='k') #neutral cells in black
                   
                ax2.plot(d[:,1], d[:,3]/d[:,2], c='r', linewidth=0.2)
                ax2.plot(d[:,1], d[:,4]/d[:,2], c='k', linewidth=0.2)
        means = np.mean(data_fracs[:,:,2:], axis=0) # first 2 cols are run number and population size, and must not be averaged
        for i in range(0,means.shape[1]-1):
            ax1.plot(data_fracs[0][:,1], means[:,i], '-', linewidth=2)
        ax1.plot(data_fracs[0][:,1], means[:,-1], '-', c='k', linewidth=2, label='neutral')
        ax1.legend(loc=1)
        ax1.set_xscale('log')
        for i in range(0,means.shape[1]-1):
            ax2.plot(data_fracs[0][:,1], means[:,i]/means[:,0], '-', linewidth=2)
        ax2.plot(data_fracs[0][:,1], means[:,-1]/means[:,0], '-', c='k', linewidth=2, label='$N_n$/$N_a$')
        ax2.legend(loc=1)
        ax2.set_xscale('log')
        
        if only_means:
            fig1.savefig(os.path.join(
                output_dir, f.split('.')[0])+'_means.png')
            fig2.savefig(os.path.join(
                output_dir, f.split('.')[0])+'_ratios_means.png')
        else:
            fig1.savefig(os.path.join(
                output_dir, f.split('.')[0])+'.png')
            fig2.savefig(os.path.join(
                output_dir, f.split('.')[0])+'_ratios.png')
        plt.close()

        # explicitly clear memory
        del data_fracs
        gc.collect()
    
    return None

if __name__ == '__main__':
        import argparse

        parser = argparse.ArgumentParser(description='Generate plots using cell fraction data')
        parser.add_argument('--input_dir', required=True,
                        help='path to simulation results')
        parser.add_argument('--output_dir', required=True,
                        help='path to output dir for saving plots')
        parser.add_argument('--runs', required=False, type=int, default=20, help='number of runs')
        parser.add_argument('--num_cells', required=False, type=int, default=100, help='number of cells in final population')
        parser.add_argument('--only_means', required=False, type=bool, default=True, help='only plot mean trajectories')
        args = parser.parse_args()
        plot_cell_fractions(input_dir=args.input_dir, output_dir=args.output_dir, runs=args.runs, num_cells=args.num_cells, only_means=args.only_means)
