
def parse_file(filename):
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
        filecontent = f.read()
        # for the last sequence processing > is added to the filecontent
        for line in filecontent.splitlines() + ['>']:
            if line.startswith('>'):
                if name:
                    joined_text = ''.join(seq_struct_text)
                    second_part = joined_text[len(joined_text) // 2:]
                    if second_part.count('.') + second_part.count('(') + \
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


def convert_to_file_data(file_data, dotbracket_structures):
    new_file_data = []
    for i in range(len(dotbracket_structures)):
        sequence = file_data[i][1]
        new_sequence = []
        structure = dotbracket_structures[i]
        b = 0
        for a in range(len(structure)):
            if structure[a] == '-':
                new_sequence.append('-')
            else:
                while sequence[b] == '-':
                    b += 1
                new_sequence.append(sequence[b])
                b += 1
        new_sequence = ''.join(new_sequence)
        new_file_data.append((
            file_data[i][0], new_sequence, structure))
    return new_file_data
