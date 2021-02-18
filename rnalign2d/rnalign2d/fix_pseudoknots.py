import argparse
import string

try:
    from .common import parse_file
except SystemError:
    from common import parse_file

OPENING = ['(', '[', '{', '<' ]
OPENING.extend(string.ascii_uppercase)
CLOSING = [')', ']', '}', '>' ]
CLOSING.extend(string.ascii_lowercase)


def structure_to_representation(dotbracket_structure):
    matching_positions = {}
    matching_pairs = {')': '(', '}': '{', ']': '[', '>': '<'}
    for letter in string.ascii_lowercase:
        matching_pairs[letter] = letter.upper()
    for position, dotbracket in enumerate(dotbracket_structure):
        if dotbracket in OPENING:
            matching_positions[position] = dotbracket
        elif dotbracket in CLOSING:
            for back_position in range(position, -1, -1):
                if back_position in matching_positions:
                    if matching_positions[back_position] == \
                            matching_pairs[dotbracket]:
                        matching_positions[back_position] = position
                        break
    for key in list(matching_positions.keys()):
        matching_positions[matching_positions[key]] = key
    return matching_positions


def representation_to_structure(dotbracket_structure, matching_positions):
    level_closures = []
    level = 0
    current_level = 0
    new_structure = []
    for index, dotbracket in enumerate(dotbracket_structure):
        if dotbracket == '.':
            new_structure.append('.')
        elif dotbracket in OPENING: # ([{ etc
            if level_closures: #remove dangling None at the end
                while level_closures[-1] == None:
                    level_closures.pop()
            if level_closures: # if there are entries in the level_closures
                if not None in level_closures \
                        and level_closures[level] > matching_positions[index]:
                    # if this is the same level of brackets - just change
                    # level_closures - to be used in case of proper CLOSING
                    level_closures[level] = matching_positions[index]
                elif None in level_closures:
                    #set level to where is None or below if possible
                    # (not nested pseudoknots on lower level)
                    level = level_closures.index(None)
                    stop = False
                    while level > 0 and not stop:
                        if level_closures[level-1] > matching_positions[index]:
                            level -= 1
                            current_level -= 1
                        else:
                            stop = True
                    # for selected level set matching_positions
                    level_closures[level] = matching_positions[index]
                else:
                    # if there is not None in the level_closures and it is
                    # pseudoknot - increase level and add
                    # proper matching_positions
                    level = len(level_closures)
                    level_closures.append(matching_positions[index])
            else:
                level_closures.append(matching_positions[index])
            # for given level add proper sign to new secondary structure
            new_structure.append(OPENING[level])
        elif dotbracket in CLOSING: #)]} etc
            if index in level_closures:
                #if it is time to close brackets
                # remove current level from the list
                current_level = level_closures.index(index)
                level_closures[current_level] = None
            # add reverse brackets:
            new_structure.append(
                CLOSING[OPENING.index(
                    new_structure[matching_positions[index]])])
    return ''.join(new_structure)


def process_file(filename_in, filename_out):
    sequences = parse_file(filename_in)
    result = []
    for name, sequence, structure in sequences:
        result.append(name)
        result.append(sequence)
        if '[' in structure:
            new_structure = representation_to_structure(
                structure, structure_to_representation(structure))
            result.append(new_structure)
        else:
            result.append(structure)
    with open(filename_out, 'w') as f_out:
        f_out.write(''.join(result))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="Input file (dot bracket)",
                        required=True)
    parser.add_argument("-o", help="Output file (dot bracket)",
                        required=True)
    args = parser.parse_args()
    process_file(args.i, args.o)


if __name__ == '__main__':
    main()
