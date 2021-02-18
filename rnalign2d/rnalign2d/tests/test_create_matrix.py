import os

import pytest

from rnalign2d.create_matrix import create_matrix


@pytest.mark.parametrize("mode, filename", [
    ("simple", "simple_matrix"),
    ("pseudo", "pseudo_matrix")])
def test_create_matrix(mode, filename):
    result = create_matrix(5, 2, -10, -8, 5, 5, mode=mode)
    file = os.path.normpath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'data', filename))
    reference = open(file, 'r').read()
    assert result == reference
