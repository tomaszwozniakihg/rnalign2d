import os

from rnalign2d.rm_mod import unmodify_file


def test_calculate_alignment_from_file():
    filename = os.path.normpath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'data', 'test_dot_bracket'))
    unmodify_file(filename, 'out_filename')
    result = open('out_filename', 'r').read()
    file = os.path.normpath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'data', 'cleared'))
    reference = open(file, 'r').read()
    os.remove('out_filename')
    assert result == reference
