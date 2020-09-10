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
      - see the documentation of the Test class below on how to set up a test case
    """
    tests = {}

    tests['tube_pressure'] = [Test('pressure', 0, 0, -1, 11005.30965, 1.0e-7, 'point'),
                              Test('pressure', 0, -1, -1, 10000.0, 1.0e-8, 'point'),
                              Test('flow', 0, -1, -1, 100.0, 1.0e-16, 'point'),
                              Test('area', 0, -1, -1, 1.0, 1.0e-5, 'point')]

    tests['tube_pressure_wave'] = [Test('pressure', 0, 0, -1, 10000.0, 1.0e-8, 'point'),
                                   Test('pressure', 0, -1, -1, 9086.52306835, 1.0e-4, 'point'),
                                   Test('flow', 0, -1, -1, 90.8652306835, 1.0e-4, 'point'),
                                   Test('area', 0, -1, -1, 1.0, 1.0e-5, 'point')]

    tests['tube_rcr'] = [Test('pressure', 0, 0, -1, 11005.30965, 1.0e-7, 'point'),
                         Test('pressure', 0, -1, -1, 10000.0, 1.0e-8, 'point'),
                         Test('flow', 0, -1, -1, 100.0, 1.0e-16, 'point'),
                         Test('area', 0, -1, -1, 1.0, 1.0e-5, 'point')]

    tests['tube_rcr_Pd'] = [Test('pressure', 0, 0, -1, 12005.30965, 1.0e-7, 'point'),
                            Test('pressure', 0, -1, -1, 11000, 1.0e-8, 'point'),
                            Test('flow', 0, -1, -1, 100.0, 1.0e-16, 'point'),
                            Test('area', 0, -1, -1, 1.0, 1.0e-5, 'point')]

    tests['tube_r'] = [Test('pressure', 0, 0, -1, 11005.30965, 1.0e-7, 'point'),
                       Test('pressure', 0, -1, -1, 10000.0, 1.0e-8, 'point'),
                       Test('flow', 0, -1, -1, 100.0, 1.0e-16, 'point'),
                       Test('area', 0, -1, -1, 1.0, 1.0e-5, 'point')]

    tests['tube_r_Pd'] = [Test('pressure', 0, 0, -1, 12005.30965, 1.0e-7, 'point'),
                          Test('pressure', 0, -1, -1, 11000.0, 1.0e-8, 'point'),
                          Test('flow', 0, -1, -1, 100.0, 1.0e-16, 'point'),
                          Test('area', 0, -1, -1, 1.0, 1.0e-5, 'point')]

    tests['tube_r_stab'] = [Test('pressure', 0, 0, -1, 11005.30965, 1.0e-7, 'point'),
                            Test('pressure', 0, -1, -1, 10000.0, 1.0e-8, 'point'),
                            Test('flow', 0, -1, -1, 100.0, 1.0e-16, 'point'),
                            Test('area', 0, -1, -1, 1.0, 1.0e-5, 'point')]

    tests['bifurcation_P'] = [Test('pressure', 0, 0, -1, 4039.45953118937, 1e-5, 'point'),
                              Test('pressure', 0, -1, -1, 4026.67220709878, 1e-5, 'point'),
                              Test('pressure', 1, 0, -1, 4026.67220709878, 1e-5, 'point'),
                              Test('pressure', 2, 0, -1, 4026.67220709878, 1e-5, 'point'),
                              Test('pressure', 1, -1, -1, 4000.00, 1e-5, 'point'),
                              Test('pressure', 2, -1, -1, 4000.00, 1e-5, 'point'),
                              Test('flow', 1, -1, -1, 3.9925, 1e-6, 'point'),
                              Test('flow', 2, -1, -1, 3.9925, 1e-6, 'point')]

    tests['bifurcation_R'] = [Test('pressure', 0, 0, -1, 3997.46433118937, 1e-5, 'point'),
                              Test('pressure', 0, -1, -1, 3984.67700709878, 1e-5, 'point'),
                              Test('pressure', 1, 0, -1, 3984.67700709878, 1e-5, 'point'),
                              Test('pressure', 2, 0, -1, 3984.67700709878, 1e-5, 'point'),
                              Test('pressure', 1, -1, -1, 3958.0048, 1e-5, 'point'),
                              Test('pressure', 2, -1, -1, 3958.0048, 1e-5, 'point'),
                              Test('flow', 1, -1, -1, 3.9925, 1e-5, 'point'),
                              Test('flow', 2, -1, -1, 3.9925, 1e-5, 'point')]

    tests['bifurcation_R_stab'] = [Test('pressure', 0, 0, -1, 3997.46433118937, 1e-6, 'point'),
                                   Test('pressure', 0, -1, -1, 3984.67700709878, 1e-6, 'point'),
                                   Test('pressure', 1, 0, -1, 3984.67700709878, 1e-6, 'point'),
                                   Test('pressure', 2, 0, -1, 3984.67700709878, 1e-6, 'point'),
                                   Test('pressure', 1, -1, -1, 3958.0048, 1e-7, 'point'),
                                   Test('pressure', 2, -1, -1, 3958.0048, 1e-7, 'point'),
                                   Test('flow', 1, -1, -1, 3.9925, 1e-6, 'point'),
                                   Test('flow', 2, -1, -1, 3.9925, 1e-6, 'point')]

    tests['bifurcation_RCR'] = [Test('pressure', 0,  0, np.arange(100, 200), 123878.022943, 1e-7, 'mean'),
                                Test('pressure', 0,  0, np.arange(100, 200), 168182.372624, 1e-7, 'max'),
                                Test('pressure', 0,  0, np.arange(100, 200), 89237.6441223, 1e-7, 'min'),
                                Test('pressure', 1, -1, np.arange(100, 200), 123855.677783, 1e-7, 'mean'),
                                Test('pressure', 1, -1, np.arange(100, 200), 171598.373528, 1e-7, 'max'),
                                Test('pressure', 1, -1, np.arange(100, 200), 87624.0897929, 1e-7, 'min'),
                                Test('pressure', 2, -1, np.arange(100, 200), 123855.677783, 1e-7, 'mean'),
                                Test('pressure', 2, -1, np.arange(100, 200), 171598.373528, 1e-7, 'max'),
                                Test('pressure', 2, -1, np.arange(100, 200), 87624.0897929, 1e-7, 'min'),
                                Test('flow', 0,  0, np.arange(100, 200), 7.557147487534, 1e-7, 'mean'),
                                Test('flow', 1, -1, np.arange(100, 200), 3.839024865141, 1e-7, 'mean'),
                                Test('flow', 1, -1, np.arange(100, 200), 24.0553490482, 1e-7, 'max'),
                                Test('flow', 1, -1, np.arange(100, 200), -3.35028925666, 1e-7, 'min')]

    tests['bifurcation_RCR_staticFunc'] = [Test('pressure', 0, 0, -1, 3997.46433118937, 1e-6, 'point'),
                                           Test('pressure', 0, -1, -1, 3984.67700709878, 1e-6, 'point'),
                                           Test('pressure', 1, 0, -1, 3984.67700709878, 1e-6, 'point'),
                                           Test('pressure', 2, 0, -1, 3984.67700709878, 1e-6, 'point'),
                                           Test('pressure', 1, -1, -1, 3958.0048, 1e-7, 'point'),
                                           Test('pressure', 2, -1, -1, 3958.0048, 1e-7, 'point'),
                                           Test('flow', 1, -1, -1, 3.9925, 1e-6, 'point'),
                                           Test('flow', 2, -1, -1, 3.9925, 1e-6, 'point')]

    return tests


class Test:
    """
    Class to define (and check) test cases
    """

    def __init__(self, field, seg, node, time, res, tol, fun):
        """
        Args:
            field: field to check ('flow', 'pressure', 'area', 'wss', 'Re')
            seg: segment to check
            node: FE nodes to check (usually 0 or -1)
            time: individual time step to check (usually -1) or list of time steps (e.g. np.arange(5,10))
            res: result to check (float)
            tol: relative tolerance for result check
            fun: type of result comparison to perform in time ('point', 'mean', 'max', 'min')
                 specifiy an interval for 'time' when chosing 'mean', 'max', or 'min'
        """
        # sanity checks
        if field not in ['flow', 'pressure', 'area', 'wss', 'Re']:
            raise ValueError('Field ' + field + ' unknown. Please select from flow, pressure, area, wss, Re')
        if fun not in ['point', 'mean', 'max', 'min']:
            raise ValueError('Function ' + fun + ' unknown. Please select from point, mean, max, min')
        if not np.isscalar(time) and fun == 'point':
            raise ValueError('Specify a single time point when selecting result type ' + fun)
        if np.isscalar(time) and (fun == 'mean' or fun == 'max' or fun == 'min'):
            raise ValueError('Specify a time interval when selecting result type ' + fun)

        self.field = field
        self.seg = seg
        self.node = node
        self.time = time
        self.res = res
        self.tol = tol
        self.fun = fun

    def check(self, results):
        """
        Perform the actual result check
        """
        # read result from svOneDSolver
        res = self.read_result(results)

        # calculate relative difference
        diff = np.abs((res - self.res) / self.res)

        # check if difference in results is larger than given tolerance
        if diff > self.tol:
            return self.print_err(res, diff)

        return False

    def read_result(self, results):
        """
        Read results and select function
        """
        # extract result
        res = results[self.field][self.seg][self.node, self.time]

        # select result type
        if self.fun == 'point':
            return res
        elif self.fun == 'mean':
            return np.mean(res)
        elif self.fun == 'max':
            return np.max(res)
        elif self.fun == 'min':
            return np.min(res)
        else:
            raise ValueError('Unknown result type ' + self.fun)

    def print_err(self, res, diff):
        """
        Create error string for user
        """
        err = 'Test failed. ' + self.field + ' in segment ' + str(self.seg)
        err += ', node ' + str(self.node) + ', time ' + str(self.time) + '. expected: ' + str(self.res)
        err += '. got: ' + str(res) + '. abs rel diff: ' + str(diff) + ' > ' + str(self.tol)
        return err


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

    # convert results to numpy array
    for f in res.keys():
        for s in res[f].keys():
            res[f][s] = np.array(res[f][s])

    return res


def run_check(results, result_checks):
    """
    Check the results of a test
    """
    # loop all results
    for test in result_checks:
        err = test.check(results)
        if err:
            # test failed
            return err

    # all tests passed
    return False


def run_test(build_dir, test_dir, name, check):
    """
    Run a test case and check the results
    """
    # name of svOneDSolver executable
    exe = os.path.join(build_dir, 'bin', 'OneDSolver')

    # name of input file
    inp = os.path.join(test_dir, name + '.in')

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


def main(solver='_skyline'):
    """
    Loop over all test cases and check if all results match
    """
    # set paths
    if 'BUILD_DIR' not in os.environ and 'TEST_DIR' not in os.environ:
        # run locally
        fpath = os.path.dirname(os.path.realpath(__file__))
        build_dir = os.path.join(fpath, '..', 'build' + solver)
        test_dir = fpath
    else:
        # run on Travis
        build_dir = os.environ['BUILD_DIR']
        test_dir = os.environ['TEST_DIR']

    # get test cases
    try:
        test_cases = get_tests()
    except Exception as err:
        print(err)
        return True

    # loop all test cases
    for name, check in test_cases.items():
        print('Running test ' + name)
        err = run_test(build_dir, test_dir, name, check)

        # check if errors occured
        if err:
            print(err)
            return True
        else:
            print('Test passed')

    # no tests failed
    return False


if __name__ == '__main__':
    # tests fail
    if main():
        exit(1)

    # tests passs
    else:
        exit(0)
