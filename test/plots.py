#!/usr/bin/env python
# coding=utf-8

import pdb
import os
import subprocess
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt

from run_tests import read_results_1d


def run_local(name):
    """
    Runs test case and returns results
    Args:
        name: name of test case
    Returns:
        dictionary res[result field][segment id][node, time step]
    """
    # name of svOneDSolver executable
    solver = '_skyline'

    # path to this folder
    fpath = os.path.dirname(os.path.realpath(__file__))

    # patth to executable
    build_dir = os.path.join(fpath, '..', 'build' + solver)
    exe = os.path.join(build_dir, 'bin', 'OneDSolver')

    # name of input file
    inp = os.path.join(fpath, name + '.in')

    # run simulation
    subprocess.check_output([exe, inp])

    # return output
    return read_results_1d('.', 'results_' + name + '_seg*')


def plot_tube_r():
    """
    Run tube_r.in and plot flow/pressure over the length of the tube
    """
    # run simulation
    res = run_local('tube_r')

    # number of finite elements
    n_fem = 50

    # length of tube
    l_tube = 10

    # axial coordinate
    x = np.linspace(0, l_tube, n_fem + 1)

    # create plot
    fig, ax1 = plt.subplots(figsize=(8, 4), dpi=300)

    # plot pressure
    ax1.plot(x, res['pressure'][0], 'b-')
    ax1.set_xlabel('Axial coordinate')
    ax1.set_ylabel('Pressure')
    ax1.grid(True)

    # plot flow
    ax2 = ax1.twinx()
    ax2.plot(x, res['flow'][0], 'r-')
    ax2.set_ylabel('Flow')

    # save figure
    fig.savefig('tube_r_results.png', bbox_inches='tight')
    plt.close(fig)


def plot_bifurcation_rcr():
    """
    Run bifurcation_rcr.in and plot flow/pressure over time at different locations
    """
    # run simulation
    res = run_local('bifurcation_RCR')

    # number of exported time steps per cardiac cycle
    n_time = 100

    # time for one cardiac cycle
    t_max = 1.087

    # time vector
    t = np.linspace(0, t_max, n_time)

    # plot locations [segment, FEM node] and their names
    locations = [[0, 0], [0, 5], [0, -1], [1, 5], [1, -1]]

    # location names
    names = ['inlet', 'mid', 'junction', 'il mid', 'il outlet']

    # fields to plot
    fields = ['flow', 'pressure']

    # create plot
    fig, ax = plt.subplots(2, 1, figsize=(8, 8), dpi=300)

    # plot results in last cardiac cycle over time at selected locations
    for i, f in enumerate(fields):
        for k in locations:
            ax[i].plot(t, res[f][k[0]][k[1], -n_time:])
        ax[i].set_xlabel('Time [s]')
        ax[i].set_ylabel(f.capitalize())
        ax[i].legend(names)
        ax[i].grid(True)

    # save figure
    fig.savefig('bifurcation_rcr_results.png', bbox_inches='tight')
    plt.close(fig)


def main():
    """
    Generate all plots
    """
    plot_tube_r()
    plot_bifurcation_rcr()


if __name__ == '__main__':
    main()
