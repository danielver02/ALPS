import numpy as np
import pytest
import os

def read_roots(path=".", filename="test_kpar_fast.scan_kpara_1.root_1"):

    path_to_file = os.path.join(path,filename)
    assert os.path.exists(path_to_file)
    data = np.loadtxt(path_to_file)
    return data

@pytest.mark.parametrize("path", ["."])
@pytest.mark.parametrize("column", range(4))
def test_check_outputs(path, column):

    ref_path = "../solution"
    
    ref_data = read_roots(ref_path)
    test_data = read_roots(path)

    np.testing.assert_almost_equal(ref_data[:, column],
                                   test_data[:, column],
                                   decimal = 4,
                                   verbose = True)
