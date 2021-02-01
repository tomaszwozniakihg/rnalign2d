import os

import pytest

from rnalign2d.rnalign2d import convert_sequence, revert_sequence, \
    remove_modifications, add_original_modifications, calculate_alignment, \
    calculate_alignment_from_file


@pytest.mark.parametrize("sequence, secondary_structure, mode, result", [
    ('AGUCCCC', '....([)', 'simple', 'DVPIGHL'),
    ('AGUCCCC', '....([)', 'pseudo', 'AAAACED'),
])
def test_convert_sequence(sequence, secondary_structure, mode, result):
    my_result = convert_sequence(sequence, secondary_structure, mode)
    assert result == my_result


@pytest.mark.parametrize("sequence, original_sequence, mode, result", [
    ('DVPIGHL', 'AGUCCCC', 'simple', ('AGUCCCC', '....([)')),
    ('AAAACED', 'AGUCCCC', 'pseudo', ('AGUCCCC', '....([)')),
    ('DVP-IGHL', 'AGUCCCC', 'simple', ('AGU-CCCC', '...-.([)')),
    ('AAA-ACED', 'AGUCCCC', 'pseudo', ('AGU-CCCC', '...-.([)')),
])
def test_revert_sequence(sequence, original_sequence, mode, result):
    my_result = revert_sequence(sequence, original_sequence, mode)
    assert result == my_result


@pytest.mark.parametrize("sequence, result", [
    ('Ab<%GK', 'AACCGG')
])
def test_remove_modifications(sequence, result):
    my_result = remove_modifications(sequence)
    assert result == my_result


@pytest.mark.parametrize("sequence, original_sequence, result", [
    ('AAC--CGG', 'Ab<%GK', 'Ab<--%GK')
])
def test_add_original_modifications(sequence, original_sequence, result):
    my_result = add_original_modifications(sequence, original_sequence)
    assert result == my_result


@pytest.mark.parametrize("sequences, mode, result", [
    ([('>tdbR00000365',
       'AAAUAUGA"GCGAUUUAUUGCAAPUAGPUUCGACCUAAUCUUAGGUGAAAUUCACCCAPAUUUUCCA',
       '(((((((..((((....)))).(((((.......)))))....((((.....)))))))))))....'),
      ('>tdbR00000030',
       'AAAAAAUU"GUUUAAUCAAAAACCPPAGUAUGUC6AACUAAAAAAAUUAGAUCAUCUAAUAPPUUU'
       'UACCA',
       '(((((((..((((......)))).(((((.......)))))....(((((.....)))))))))))'
       ')....')],
     'simple',
     [('>tdbR00000365',
       'AAAUAUGA"GCGA--UUUAUUGCAAPUAGPUUCGACCUAAUCUUAGGUG--AAAUUCACCCAPAUU'
       'UUCCA',
       '(((((((..((((--....)))).(((((.......)))))....((((--.....))))))))))'
       ')....'),
      ('>tdbR00000030',
       'AAAAAAUU"GUUUAAUCAAAAACCPPAGUAUGUC6AACUAAAAAAAUUAGAUCAUCUAAUAPPUUU'
       'UACCA',
       '(((((((..((((......)))).(((((.......)))))....(((((.....)))))))))))'
       ')....'),]),

])
def test_calculate_alignment(sequences, mode, result):
    simple_matrix = os.path.normpath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), '..', 'data',
        'simple_matrix'))
    my_result = calculate_alignment(
        sequences, mode, simple_matrix, -12, -1)
    assert result == my_result


@pytest.mark.parametrize(
    "testfile, resultfile", [
    ('test_dot_bracket', 'reference'),
    ('test_dot_bracket_multiline', 'reference'),
])
def test_calculate_alignment_from_file(testfile, resultfile):
    filename = os.path.normpath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'data', testfile))
    mode = 'simple'
    matrix = os.path.normpath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), '..', 'data',
        'simple_matrix'))
    gapopen = -12
    gapextend = -1
    calculate_alignment_from_file(
            filename, 'out_filename', mode, matrix, gapopen, gapextend)
    result = open('out_filename', 'r').read()
    file = os.path.normpath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'data', resultfile))
    reference = open(file, 'r').read()
    os.remove('out_filename')
    assert result == reference
