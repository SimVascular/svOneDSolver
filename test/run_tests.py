#!/usr/bin/env python
# coding=utf-8

import csv
import glob
import os
import pdb
import re
import subprocess
import numpy as np
from collections import defaultdict


def get_tests():
    """
    Add new test cases here

    How to define a new test case <NAME>:
      - create new file with name <NAME>.in with the following naming conventions:
        MODEL results_<NAME>_
        SEGMENT seg0 ...
        SEGMENT seg1 ...
      - add test case to the dictionary below
      - all dictionary values are lists, so add as many results as you like

    tests[<NAME>]['field'] = [fields to check ('flow', 'pressure', 'area', 'wss', 'Re')]
    tests[<NAME>]['seg'] = [segments to check]
    tests[<NAME>]['node'] = [FE nodes to check (usually 0 or -1)]
    tests[<NAME>]['time'] = [time steps to check (usually -1)]
    tests[<NAME>]['res'] = [results to check]
    tests[<NAME>]['tol'] = [tolerances for result check]

    """
    tests = defaultdict(dict)

    tests['tube_pressure']['field'] = ['pressure', 'pressure']
    tests['tube_pressure']['seg'] = [0, 0]
    tests['tube_pressure']['node'] = [0, -1]
    tests['tube_pressure']['time'] = [-1, -1]
    tests['tube_pressure']['res'] = [11005.30965, 10000]
    tests['tube_pressure']['tol'] = [1.0e-7, 1.0e-9]

    tests['tube_rcr']['field'] = ['pressure', 'pressure']
    tests['tube_rcr']['seg'] = [0, 0]
    tests['tube_rcr']['node'] = [0, -1]
    tests['tube_rcr']['time'] = [-1, -1]
    tests['tube_rcr']['res'] = [11005.30965, 10000]
    tests['tube_rcr']['tol'] = [1.0e-7, 1.0e-8]

    tests['tube_r']['field'] = ['pressure', 'pressure']
    tests['tube_r']['seg'] = [0, 0]
    tests['tube_r']['node'] = [0, -1]
    tests['tube_r']['time'] = [-1, -1]
    tests['tube_r']['res'] = [11005.30965, 10000]
    tests['tube_r']['tol'] = [1.0e-7, 1.0e-8]
    return tests


def read_results_1d(res_dir, name):
    """
    Read results from oneDSolver and store in dictionary
    Args:
        res_dir: directory containing 1D results
        name: pattern to search for in result file names
    Returns:
        dictionary res[result field][segment id][node, time step]
    """
    # read from files, store in dict, and remove files
    res = defaultdict(lambda: defaultdict(list))
    for field in ['flow', 'pressure', 'area', 'wss', 'Re']:
        for f_res in glob.glob(os.path.join(res_dir, name + '_' + field + '.dat')):
            with open(f_res) as f:
                for line in csv.reader(f, delimiter=' '):
                    res[field][int(re.findall(r'\d+', f_res)[-1])] += [[float(m) for m in line if m][1:]]
            os.remove(f_res)

    # convert to numpy array
    for f in res.keys():
        for s in res[f].keys():
            res[f][s] = np.array(res[f][s])

    return res


def run_check(results, c):
    """
    Check the results of a test
    """
    for field, seg, node, time, ref, tol in zip(c['field'], c['seg'], c['node'], c['time'], c['res'], c['tol']):
        res = results[field][seg][node, time]
        diff = np.abs(res - ref) / ref

        if diff > tol:
            err = 'Test failed. ' + field + ' in segment ' + str(seg) + ', node ' + str(node) + ', time ' + str(time)
            err += '. expected: ' + str(ref) + '. got: ' + str(res) + '. abs rel diff: ' + str(diff) + ' > ' + str(tol)
            return err
    else:
        return False


def run_test(name, check, solver='superlu'):
    """
    Run a test case and check the results
    """
    # name of svOneDSolver executable
    exe = os.path.join('..', 'build_' + solver, 'bin', 'OneDSolver')

    # name of input file
    inp = name + '.in'

    # run simulation
    try:
        subprocess.check_output([exe, inp])
    except subprocess.CalledProcessError as err:
        return 'Test failed. svOneDSolver returned error:\n' + err.output.decode("utf-8")

    # extract results
    try:
        res = read_results_1d('.', 'results_' + name + '_seg*')
    except Exception as err:
        return 'Test failed. Result extraction failed:\n' + str(err)

    # compare to stored results
    try:
        res = run_check(res, check)
    except Exception as err:
        return 'Test failed. Result check failed:\n' + str(err)

    return res


def main(solver='superlu'):
    """
    Loop over all test cases and check if all results match
    """
    for name, check in get_tests().items():
        print('Running test ' + name)
        err = run_test(name, check, solver)
        if err:
            print(err)
            return True
        else:
            print('Test passed')
    else:
        return False


if __name__ == '__main__':
    if main():
        exit(-1)
    else:
        exit(0)
