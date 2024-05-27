import numpy as np
import pytest
import os

def read_roots(path=".", filename="test_kpar_fast.scan_kpara_1.root_1"):

    path_to_file = os.path.join(path,filename)
    assert os.path.exists(path_to_file)
    data = np.loadtxt(path_to_file)
    return data

@pytest.mark.parametrize("filename", ["test_kpar_fast.scan_kpara_1.root_1",
                                      "test_kpar_fast.eigen_kpara_1.root_1",
                                      "test_kpar_fast.heat_kpara_1.root_1"])
def test_check_outputs(filename):
    '''
    tests that the contents of the file solution/<filename> matches the file
    ./<filename> to a relative tolerance of 1e-3. Will do this for each file
    in the parameter list.
    '''

    path = "../solution"
    ref_path = "."

    ref_data = read_roots(ref_path, filename)
    test_data = read_roots(path, filename)

    np.testing.assert_allclose(ref_data,
                               test_data,
                               rtol = 1.e-3,
                               atol = 1.e-4,
                               verbose = True)
