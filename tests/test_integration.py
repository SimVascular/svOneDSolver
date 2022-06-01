import os

import shutil
import numpy as np
import csv
import glob
import os
import re
import subprocess

RTOL_PRES = 1.0e-7
RTOL_FLOW = 1.0e-8


def run_test_case_by_name(name, testdir):
    """Run a test case by its case name.

    Args:
        name: Name of the test case.
        testdir: Directory for performing the simulation.
    """
    # file where executable has been built
    build_dir = 'build_skyline'

    # test file
    testfile = os.path.join(os.path.dirname(__file__), 'cases', name + '.in')
    testfile_tmp = os.path.join(testdir, name + '.in')

    # copy input file to temporary test dir
    shutil.copyfile(testfile, testfile_tmp)

    # run 1D simulation
    run_oned(build_dir, testfile_tmp)

    # extract results
    return read_results_1d('.', 'results_' + name + '_seg*')


def run_oned(build_dir, name):
    """
    Executes svOneDSolver with given input file
    """
    # name of svOneDSolver executable
    exe = os.path.join(build_dir, 'bin', 'OneDSolver')

    # name of input file
    inp = os.path.join('cases', name + '.in')

    # run simulation
    try:
        subprocess.check_output([exe, inp])
    except subprocess.CalledProcessError as err:
        raise RuntimeError('Test failed. svOneDSolver returned error:\n' + err.output.decode("utf-8"))


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


def get_result(results, field, seg, node, time, fun):
    """
    Read results and select function
    """
    # extract result
    res = results[field][seg][node, time]

    # select result type
    if fun == 'point':
        return res
    elif fun == 'mean':
        return np.mean(res)
    elif fun == 'max':
        return np.max(res)
    elif fun == 'min':
        return np.min(res)
    else:
        raise ValueError('Unknown result type ' + fun)


def test_tube_pressure(tmpdir):
    results = run_test_case_by_name('tube_pressure', tmpdir)
    assert np.isclose(get_result(results, 'pressure', 0, 0, -1, 'point'), 11005.30965, rtol=1.0e-7)
    assert np.isclose(get_result(results, 'pressure', 0, 0, -1, 'point'), 10000.0, rtol=1.0e-8)
    assert np.isclose(get_result(results, 'flow', 0, 0, -1, 'point'), 100.0, rtol=1.0e-16)
    assert np.isclose(get_result(results, 'area', 0, 0, -1, 'point'), 1.0, rtol=1.0e-5)

    # assert np.isclose(get_result(results, '', 0, 0, -1, 'point'), , rtol=)
