import argparse
import os
from uuid import uuid4
try:
    from .conversion import SIMPLE_CONVERSION, PSEUDOKNOT_CONVERSION
except SystemError:
    from rnalign2d.conversion import SIMPLE_CONVERSION, PSEUDOKNOT_CONVERSION


MODIFICATIONS = {
            'A': ['A', 'H', '\"', '/', '+', '*', '=', '6', 'E', '[', ':', 'I',
                  'O', '^', '`', 'b', '≠', 'ÿ', '«'],
            'C': ['C', '<', '%', 'B', 'M', '?', "'", '}', '>', '°', 'a', '¿',
                  '’'],
            'G': ['G', ';', 'K', 'L', '#', 'R', '|', '7', '(', 'Q', '8', '9',
                  'Y', 'W', 'd', '⊄'],
            'U': ['U', 'N', '{', '2', 'J', '4', '&', '1', 'S', '3', 'V', '5',
                  '!', '$', 'X', ',', ')', '~', 'D', 'P', ']', 'Z', 'T', 'F',
                  '\\', 't', '@', 'c', 'Ê', '∃', 'υ', 'Ð',
                  '_', '-', '.']}


def convert_sequence(sequence, secondary_structure, mode='simple'):
    new_sequence = []
    if mode == 'simple':
        for i, letter in enumerate(sequence):
            ss_letter = secondary_structure[i]
            if not ss_letter in '([.])':
                ss_letter = '.'
            new_sequence.append(SIMPLE_CONVERSION[letter][ss_letter])
    elif mode == 'pseudo':
        for letter in secondary_structure:
            if not letter in PSEUDOKNOT_CONVERSION.keys():
                letter = '.'
            new_sequence.append(PSEUDOKNOT_CONVERSION[letter])
    return "".join(new_sequence)


def revert_sequence(sequence, original_sequence, mode):
    new_sequence = []
    new_secondary_structure = []
    original_sequence = [y for y in original_sequence]
    original_sequence.reverse()
    for seq_letter in sequence:
        if seq_letter == '-':
            new_sequence.append('-')
            new_secondary_structure.append('-')
        else:
            if mode == 'simple':
                for letter in SIMPLE_CONVERSION:
                    for dot_bracket in SIMPLE_CONVERSION[letter]:
                        if SIMPLE_CONVERSION[letter][dot_bracket] == seq_letter:
                            new_sequence.append(letter)
                            new_secondary_structure.append(dot_bracket)
            elif mode == 'pseudo':
                new_sequence.append(original_sequence.pop())
                for dot_bracket in PSEUDOKNOT_CONVERSION:
                    if PSEUDOKNOT_CONVERSION[dot_bracket] == seq_letter:
                        new_secondary_structure.append(dot_bracket)
    return "".join(new_sequence), "".join(new_secondary_structure)


def remove_modifications(sequence):
    # use U if not found
    new_sequence = []
    for letter in sequence:
        found = False
        for key in MODIFICATIONS:
            if letter in MODIFICATIONS[key]:
                new_sequence.append(key)
                found = True
        if not found:
            print('Warning - not found:', letter, 'using U instead')
            new_sequence.append('U')
    return "".join(new_sequence)


def add_original_modifications(sequence, original_sequence):
    original_sequence = [y for y in original_sequence]
    original_sequence.reverse()
    new_sequence = []

    for letter in sequence:
        if letter == '-':
            new_sequence.append('-')
        else:
            new_sequence.append(original_sequence.pop())
    return "".join(new_sequence)


def calculate_alignment(sequences, mode, matrix, gapopen, gapextend):
    """
    1 - remove modifications
    2 - convert sequence
    3 - muscle - msa
    4 - revert sequences
    5 - add original modifications
    """
    hash = uuid4().hex
    new_file_lines = []
    for i, element in enumerate(sequences):
        name, sequence, structure = element
        unmodified_sequence = remove_modifications(sequence)
        converted_sequence = convert_sequence(
            unmodified_sequence, structure, mode)
        new_file_lines.append('>{}'.format(str(i)))
        new_file_lines.append('{}'.format(converted_sequence))
    new_file_content = "\n".join(new_file_lines)
    temp_name_in = os.path.join(
        os.getcwd(), 'temp_1_{}'.format(hash))
    temp_name_out = os.path.join(
        os.getcwd(), 'temp_2_{}'.format(hash))
    with open(temp_name_in, 'w') as f:
        f.write(new_file_content)

    command = 'muscle -in {} -out {} -matrix {} -gapopen {} ' \
              '-gapextend {} -center 0.0'.format(
        temp_name_in, temp_name_out, matrix, gapopen, gapextend)
    os.system(command)
    new_sequences = []
    with open(temp_name_out, 'r') as f:
        counter = 0
        name = None
        sequence = ''
        for line in f.readlines():
            if line.startswith('>'):
                if counter != 0 and len(line.strip()) > 0:
                    my_id = int(name.replace('>', ''))
                    original_sequence = sequences[my_id][1]
                    original_name = sequences[my_id][0]
                    new_sequence, new_structure = revert_sequence(
                        sequence, original_sequence, mode)
                    new_sequence = add_original_modifications(
                        new_sequence, original_sequence)
                    new_sequences.append(
                        (original_name, new_sequence, new_structure))
                    sequence = ''
                name = line.strip()
            else:
                sequence += line.strip()
            counter += 1
        my_id = int(name.replace('>', ''))
        original_sequence = sequences[my_id][1]
        original_name = sequences[my_id][0]
        new_sequence, new_structure = revert_sequence(
            sequence, original_sequence, mode)
        new_sequence = add_original_modifications(
            new_sequence, original_sequence)
        new_sequences.append((original_name, new_sequence, new_structure))
    os.remove(temp_name_in)
    os.remove(temp_name_out)
    return new_sequences


def calculate_alignment_from_file(
        filename, out_filename, mode, matrix, gapopen, gapextend):
    sequences = []
    name = None
    sequence = None
    structure = None
    seq_struct_text = []
    counter = 0

    vienna = False
    try:
        import RNA
        vienna = True
    except:
        pass

    with open(filename, 'r') as f:
        filecontent = f.read()
        # for the last sequence processing > is added to the filecontent
        for line in filecontent.splitlines() + ['>']:
            if line.startswith('>'):
                if name:
                    joined_text = ''.join(seq_struct_text)
                    second_part = joined_text[len(joined_text)//2:]
                    if second_part.count('.') + second_part.count('(') + \
                            second_part.count(')') + second_part.count('[') + \
                            second_part.count(']') > len(second_part) / 2:
                        sequence = joined_text[:len(joined_text)//2]
                        structure = second_part
                    else:
                        sequence = joined_text
                        if vienna:
                            fc = RNA.fold_compound(sequence)
                            structure, mfe = fc.mfe()
                        else:
                            structure = '.' * len(sequence)
                    sequences.append((name, sequence, structure))
                    seq_struct_text = []
                name = line
            else:
                seq_struct_text.append(line)
    result = calculate_alignment(sequences, mode, matrix, gapopen, gapextend)
    with open(out_filename, 'w') as f:
        for element in result:
            f.write("{}\n{}\n{}\n".format(*element))


def main():
    parser = argparse.ArgumentParser()
    simple_matrix = os.path.normpath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'data', 'simple_matrix'))
    parser.add_argument("-i", help="Input file (dot bracket)", required=True)
    parser.add_argument("-o", help="Output file (dot bracket)", required=True)
    parser.add_argument(
        "-matrix", help="Matrix for alignment", default=simple_matrix)
    parser.add_argument(
        "-mode", help="Mode for matrix creation - it can be 'simple' for "
                      "comparision with maximum of one level of pseudoknots "
                      "or 'pseudo' for multiple level of pseudoknots",
        choices=['simple', 'pseudo'], default='simple')
    parser.add_argument("-gapopen", type=int, default=-12)
    parser.add_argument("-gapextend", type=int, default=-1)
    args = parser.parse_args()
    calculate_alignment_from_file(
        args.i, args.o, mode=args.mode, matrix=args.matrix,
        gapopen=args.gapopen, gapextend=args.gapextend)


if __name__ == '__main__':
    main()
