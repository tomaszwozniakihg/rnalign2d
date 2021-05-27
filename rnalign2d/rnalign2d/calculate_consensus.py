import argparse
from collections import defaultdict
import string


def calculate_consensus(filename, raw):
    file_data = parse_file(filename, raw)
    if len(file_data) == 0:
        return 0
    dotbracket_structures = [x[2] for x in file_data]
    representations = [
        structure_to_representation(structure)
        for structure in dotbracket_structures]
    # assuming that )] [( is greater than . and it is greater than -
    minimum_count = len(file_data)// 2 + len(file_data) % 2

    consensus = []
    for position in range(len(file_data[0][2])):
        counter = defaultdict(int)
        for structure in dotbracket_structures:
            counter[structure[position]] += 1
        # calculate consensus
        found_any = False
        for key in ['(', ')', '[', ']', '{', '}', '.', '-']:
            if counter[key] >= minimum_count:
                if key not in ['.', '-']:
                    counter2 = defaultdict(int)
                    for representation in representations:
                        if position in representation:
                            counter2[representation[position]] += 1
                    can_be_db = False
                    for my_key in counter2.keys():
                        if counter2[my_key] >= minimum_count:
                            can_be_db = True
                            break
                    if can_be_db:
                        consensus.append(key)
                        found_any = True
                    else:
                        consensus.append('.')
                        found_any = True
                else:
                    consensus.append(key)
                    found_any = True
            if found_any:
                break
        if not found_any:
            consensus.append('.')

    return "".join(consensus)


def structure_to_representation(dotbracket_structure):
    matching_positions = {}
    matching_pairs = {')': '(', '}': '{', ']': '[', '>': '<'}
    for letter in string.ascii_lowercase:
        matching_pairs[letter] = letter.upper()
    for position, dotbracket in enumerate(dotbracket_structure):
        if dotbracket in ('(', '<', '[', '{') or dotbracket.isupper():
            matching_positions[position] = dotbracket
        if dotbracket in (')', '>', ']', '}') or dotbracket \
                in string.ascii_lowercase:
            for back_position in range(position, -1, -1):
                if back_position in matching_positions:
                    if matching_positions[back_position] == \
                            matching_pairs[dotbracket]:
                        matching_positions[back_position] = position
                        break
    for key in list(matching_positions.keys()):
        matching_positions[matching_positions[key]] = key
    return matching_positions


def parse_file(filename, raw):
    sequences = []
    name = None
    seq_struct_text = []

    vienna = False
    try:
        import RNA
        vienna = True
    except:
        pass

    with open(filename, 'r') as f:
        if raw:
            for line in f:
                sequences.append(('', '', line.strip()))
            return sequences
        filecontent = f.read()
        # for the last sequence processing > is added to the filecontent
        for line in filecontent.splitlines() + ['>']:
            if line.startswith('>'):
                if name:
                    joined_text = ''.join(seq_struct_text)
                    second_part = joined_text[len(joined_text) // 2:]
                    if second_part.count('.') + second_part.count('-') \
                            + second_part.count('(') + \
                            second_part.count(')') + second_part.count('[') + \
                            second_part.count(']') > len(second_part) / 2:
                        sequence = joined_text[:len(joined_text) // 2]
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
    return sequences


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='Input')
    parser.add_argument('-r', help='Raw structures', action='store_true')
    args = parser.parse_args()
    print(calculate_consensus(args.i, args.r))
