import argparse
from collections import defaultdict
import string
try:
    from .common import parse_file, convert_to_file_data
except SystemError:
    from rnalign2d.common import parse_file, convert_to_file_data


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


def structure_conservation(dotbracket_structures):
    if len(dotbracket_structures) == 0:
        return 0
    conservance = []
    for i in range(len(dotbracket_structures[0])):
        letters = defaultdict(int)
        for structure in dotbracket_structures:
            letters[structure[i]] += 1
        score = 0
        for letter in letters.keys():
            score += letters[letter] ** 2
        conservance.append(score)
    return conservance


def score_by_conservation(dotbracket_structures):
    conservation = structure_conservation(dotbracket_structures)
    return sum(conservation)/len(conservation)


def find_structural_blocks(dotbracket_structure, representation):
    blocks = []
    dotbracket_structure_len = len(dotbracket_structure)
    current_block_start = None
    for i in range(dotbracket_structure_len):
        if dotbracket_structure[i] not in ('.', '-'):
            if current_block_start == None:
                current_block_start = i
                if i + 1 < dotbracket_structure_len:
                    if dotbracket_structure[i] != dotbracket_structure[i+1]:
                        blocks.append((current_block_start, i))
                        current_block_start = None
                else:
                    blocks.append((i, i))

            else:
                if i + 1 < dotbracket_structure_len:
                    if dotbracket_structure[i] == dotbracket_structure[i+1] \
                            and dotbracket_structure[representation[i]] == \
                            dotbracket_structure[representation[i+1]]:
                        # if end of structure:
                        if dotbracket_structure_len == i + 1:
                            blocks.append((current_block_start, i+1))

                    else:
                        # if there is mismatch
                        if i+2 < dotbracket_structure_len:
                            if dotbracket_structure[i] == \
                                    dotbracket_structure[i+2] and \
                                    dotbracket_structure[representation[i]] \
                                    == dotbracket_structure[
                                        representation[i+2]]:
                                continue
                        # if there is no mismatch
                        blocks.append((current_block_start, i))
                        current_block_start = None
                else:
                    blocks.append((current_block_start, i))
    return blocks


def calculate_unusual_positions_places(
        representations, structural_blocks, max_diff=2):
    """
    returns list of unusual positions
    """
    unusual_positions_places = []
    unusual_positions_and_length = []
    for r1 in range(len(representations)):
        representation1 = representations[r1]
        blocks1 = structural_blocks[r1]
        for r2 in range(len(representations)):
            representation2 = representations[r2]
            blocks2 = structural_blocks[r2]
            for i in representation1.keys():
                if i in representation2.keys():
                    if representation1[i] != representation2[i] \
                            and abs(representation1[i] - representation2[i]) \
                                    <= max_diff:
                        min_position = min(
                            i, representation1[i],
                            representation2[i],
                            representation2[representation2[i]])
                        unusual_positions_places.append(min_position)
                        unusual_positions_and_length.append(
                            (min_position,
                             abs(representation1[i] - representation2[i])))
            for block1 in blocks1:
                for a in range(*block1):
                    for block2 in blocks2:
                        if a in range(*block2):
                            if len(range(*block1)) == len(range(*block2)):
                                if block1[0] != block2[0]:
                                    data = range(block1[0], block1[1] + 1) \
                                        if block1[0] < block2[0] else \
                                        range(block2[0], block2[1] + 1)
                                    diff = abs(block1[0] - block2[0])
                                    for position in data:
                                        unusual_positions_places.append(
                                            position)
                                        unusual_positions_and_length.append(
                                            (position, diff))
    return sorted(list(set(unusual_positions_places))), \
           sorted(list(set(unusual_positions_and_length)))


def move_1_2nt_gaps(
        dotbracket_structures, offset=0, max_diff=5, multi_score=1.01,
        offset_time=0):
    representations = [
        structure_to_representation(structure)
        for structure in dotbracket_structures]
    structural_blocks = [
        find_structural_blocks(dotbracket_structure, representation)
        for dotbracket_structure, representation in
        zip(dotbracket_structures, representations)]
    unusual_positions_places, unusual_positions_and_length = \
        calculate_unusual_positions_places(
            representations, structural_blocks, max_diff)
    unusual_and_blocks_status = []
    for unusual_position, how_many_nt in unusual_positions_and_length:
        if unusual_position < offset:
            continue
        #for each structure find out if '-' is inside, pre or post
        # given block - check also regions between blocks
        for structure_no, single_structure_blocks in \
                enumerate(structural_blocks):
            for block_no, block in enumerate(single_structure_blocks):
                gap_in = False
                if '-' in dotbracket_structures[structure_no][
                          block[0]:block[1]]:
                    gap_in = True
                left = 0

                right = len(dotbracket_structures[structure_no]) - 1
                if block_no < len(single_structure_blocks) - 1:
                    right = single_structure_blocks[block_no+1][0] - 1
                if unusual_position in range(block[0]+1, block[1]):
                    unusual_and_blocks_status.append(
                        [block_no, ('right', gap_in, left, right)])
                    break
                elif unusual_position == block[0]:
                    unusual_and_blocks_status.append(
                        [block_no, ('left', gap_in, left, right)])
                    break
        # calculate consensus
        consensus_dict = defaultdict(int)
        # in rare case there is no consensus
        if unusual_and_blocks_status == []:
            continue
        for u_b_stat in unusual_and_blocks_status:
            consensus_dict[u_b_stat[1]] += 1
        max_value = 0
        consensus = None
        for key in consensus_dict:
            if consensus_dict[key] > max_value:
                max_value = consensus_dict[key]
                consensus = key
        left_or_right = consensus[0]

        score_pre = score_by_conservation(dotbracket_structures)

        offset_time += 1
        if offset_time <= 2:
            solution = fix_one_place_constant_dist(
                dotbracket_structures=dotbracket_structures,
                structural_blocks=structural_blocks,
                position=unusual_position,
                unusual_positions_places=unusual_positions_places,
                representations=representations,
                how_many_nt=how_many_nt)
            if solution:
                return move_1_2nt_gaps(
                    solution, offset=offset, offset_time=offset_time)
        solution = fix_one_place(
            dotbracket_structures=dotbracket_structures,
            position=unusual_position,
            left_or_right=left_or_right,
            unusual_positions_places=unusual_positions_places,
            representations=representations,
            how_many_nt=how_many_nt)
        if not solution:
            return move_1_2nt_gaps(dotbracket_structures, offset=offset)
        score_post = score_by_conservation(solution)

        offset = consensus[3]
        if score_post * multi_score >= score_pre:
            return move_1_2nt_gaps(solution, offset=offset)
        else:
            return move_1_2nt_gaps(dotbracket_structures, offset=offset)
    # if no further change is possible
    return dotbracket_structures


def find_counter_start_end(dotbracket_structures, position,
        unusual_positions_places, representations):
    # find end position good
    end_position = position
    while True:
        if end_position + 1 in unusual_positions_places:
            end_position += 1
        else:
            break

    counter_start = 1000000000
    counter_end = 1000000000
    for dotbracket, representation in \
            zip(dotbracket_structures, representations):
        if dotbracket[position] not in ('.', '-'):
            counter_start = min(counter_start, representation[position])
            my_position = end_position
            while dotbracket[my_position] in ('.', '-'):
                my_position -= 1
            counter_end = min(counter_end, representation[my_position])
        if dotbracket[end_position] not in ('.', '-'):
            counter_end = min(counter_end, representation[end_position])
            my_position = position
            while dotbracket[my_position] in ('.', '-'):
                my_position += 1
            counter_start = min(counter_start, representation[my_position])
    return counter_end, counter_start


def fix_one_place_constant_dist(
        dotbracket_structures, structural_blocks, position,
        unusual_positions_places, representations, how_many_nt):
    def _fix_constant_len_blocks(x_position):
        block_len = 0
        block_start = 100000000
        block_max_start = 0
        same_len = True
        block_ids = []
        for blocks in structural_blocks:
            for block_id, block in enumerate(blocks):
                if x_position + how_many_nt in range(block[0], block[1] + 1):
                    block_ids.append(block_id)
                    block_start = min(block[0], block_start)
                    block_max_start = max(block[0], block_max_start)
                    if block_len == 0:
                        block_len = block[1] - block[0]
                    else:
                        if block_len != block[1] - block[0]:
                            same_len = False
        if same_len and block_start != block_max_start \
                and len(structural_blocks) == len(block_ids):
            # if block of the same len and not starting in the same position
            # and overlapping position (just to simplify)
            shift_table = [block_max_start -
                           structural_blocks[i][block_ids[i]][0]
                           for i in range(len(structural_blocks))]

            # mock
            how_many_nt_by_structure = shift_table
            how_many_nt_max = max(shift_table)
            new_structures = move_structures(
                dotbracket_structures, start_position, end_position,
                'left', how_many_nt_max, how_many_nt_by_structure)
            return new_structures

    #calculate position
    end_position = position
    start_position = position
    while True:
        if end_position + 1 in unusual_positions_places:
            end_position += 1
        else:
            break

    # calculate score for the positon and counter_position to determine what
    # to move
    # change position to counter position if necessary
    counter_start = 0
    counter_end = 10000000
    for dotbracket, representation in \
            zip(dotbracket_structures, representations):
        if dotbracket[position] not in ('.', '-') \
                and dotbracket[end_position] not in ('.', '-'):
            counter_start = max(representation[position], counter_start)
            counter_end = min(representation[end_position], counter_end)

    x_position = start_position
    result = _fix_constant_len_blocks(x_position)
    if result:
        return result

    x_position = counter_end
    result = _fix_constant_len_blocks(x_position)
    if result:
        return result
    return None


def fix_one_place(dotbracket_structures, position, left_or_right,
                  unusual_positions_places, representations, how_many_nt):
    """
    left: 
    -((-(((-
    .-(((((.
    right:
    -((-(((-
    .(((((-.
    """

    def _get_slice_of_dotbracket(dotbracket_structures, start, end):
        """
        :param dotbracket_structures: list of dot-bracket secondary structures
        :param start: starting point - where to cut
        :param end: ending point - where to cut
        :return: list of slices from dot-bracket secondary structures
        """
        return [x[start:end] for x in dotbracket_structures]

    # find end position good
    end_position = position
    start_position = position
    while True:
        if end_position + 1 in unusual_positions_places:
            end_position += 1
        else:
            break

    counter_end, counter_start = find_counter_start_end(
        dotbracket_structures, position, unusual_positions_places,
        representations)

    dotbracket_structures_orig = _get_slice_of_dotbracket(
        dotbracket_structures, position, end_position+1)
    dotbracket_structures_counter = _get_slice_of_dotbracket(
        dotbracket_structures, counter_end, counter_start+1)
    score_1 = score_by_conservation(dotbracket_structures_orig)
    score_2 = score_by_conservation(dotbracket_structures_counter)

    # check if change should be in the counter place
    if score_1 > score_2:
        # prepare for reverse
        start_position = counter_end
        end_position = counter_start
        left_or_right = 'left' if left_or_right == 'right' else 'right'

    # create two, or more groups
    groups = defaultdict(list)
    detailed_groups = defaultdict(list)
    for dotbracket, representation, structure_no in zip(
            dotbracket_structures, representations,
            range(len(dotbracket_structures))):
        if position in representation:
            groups[representation[position]].append(structure_no)
            detailed_groups[representation[position]].append(structure_no)
        else:
            groups['position'].append(structure_no)
            found = False
            for i in range(1, how_many_nt + 1):
                if position + i in representation:
                    detailed_groups[representation[position + i] + i].append(
                        structure_no)
                    found = True
                    break
                if position - i in representation:
                    detailed_groups[representation[position - i] - i].append(
                        structure_no)
                    found = True
                    break
            if not found:
                detailed_groups['position'].append(structure_no)

    #check group, that is more to the left
    groups = detailed_groups
    names = list(groups.keys())
    if 'position' in names:
        names.remove('position')

    sorted_group_keys = sorted(
        names, reverse=True if left_or_right=='left' else False)
    # if 'position' - then add it to the first group
    if 'position' in groups.keys():
        groups[sorted_group_keys[0]].extend(groups['position'])
        del groups['position']

    # if distance is acceptable - then do it - in case of lack of
    # secondary structure element - distance may be huge
    distances = [
        abs(sorted_group_keys[0] - x)
        if abs(sorted_group_keys[0] - x) <= how_many_nt
        else 0 for x in sorted_group_keys]
    how_many_nt_by_structure = []
    for i in range(len(dotbracket_structures)):
        for key_no, key in enumerate(sorted_group_keys):
            if i in groups[key]:
                how_many_nt_by_structure.append(distances[key_no])
                break


    how_many_nt_max = abs(sorted_group_keys[0] - sorted_group_keys[-1])
    new_structures = move_structures(
        dotbracket_structures, start_position, end_position, left_or_right,
        how_many_nt_max, how_many_nt_by_structure)
    return new_structures


def move_structures(dotbracket_structures, start_position, end_position,
                    left_or_right, how_many_nt, how_many_nt_by_structure):
    able_to_move = True

    def _find_right_gap_indexes(
            dotbracket_structures, start_seek_position, end_seek_position,
            seek_step, how_many_nt_by_structure):
        able_to_move = True
        right_gap_indexes = defaultdict(list)
        for structure_id in range(len(dotbracket_structures)):
            if not able_to_move:
                break
            how_may_to_find = how_many_nt_by_structure[structure_id]
            if not how_may_to_find:
                continue
            structure = dotbracket_structures[structure_id]
            for position in range(
                    start_seek_position, end_seek_position, seek_step):
                letter = structure[position]
                if letter == '-':
                    right_gap_indexes[structure_id].append(position)
                    how_may_to_find -= 1
                    if not how_may_to_find:
                        break
                elif letter != '.':
                    able_to_move = False
                    break
        return right_gap_indexes, able_to_move

    new_structures = []
    if left_or_right == 'left':
        right_gap_indexes, able_to_move = _find_right_gap_indexes(
            dotbracket_structures, end_position+1,
            len(dotbracket_structures[0]), 1, how_many_nt_by_structure)
        if able_to_move:
            new_structures = []
            for structure_index, structure in enumerate(dotbracket_structures):
                how_many_nt_structure = how_many_nt_by_structure[
                    structure_index]
                new_structure = structure[:start_position] \
                                + '-' * how_many_nt_structure

                index_start = start_position
                for index in right_gap_indexes[structure_index]:
                    index_end = index
                    new_structure += structure[index_start:index_end]
                    index_start = index_end + 1
                new_structure += structure[index_start:]
                new_structures.append(new_structure)
            return new_structures
    else:
        right_gap_indexes, able_to_move = _find_right_gap_indexes(
            dotbracket_structures,
            start_position+how_many_nt-1, 0, -1, how_many_nt_by_structure)
        if able_to_move:
            for structure_index, structure in enumerate(dotbracket_structures):
                how_many_nt_structure = how_many_nt_by_structure[
                    structure_index]
                new_structure = ''
                index_start = 0
                for index in sorted(right_gap_indexes[structure_index]):
                    index_end = index
                    new_structure += structure[index_start:index_end]
                    index_start = index_end + 1
                new_structure += \
                    structure[index_start:end_position+1+how_many_nt_structure]\
                    + '-' * how_many_nt_structure \
                    + structure[end_position+1+how_many_nt_structure:]
                new_structures.append(new_structure)
            return new_structures

    # if not possible to move - just add - in both strands
    # in respective positions
    for structure_index, structure in enumerate(dotbracket_structures):
        how_many_nt_structure = how_many_nt_by_structure[structure_index]
        residual_nt = how_many_nt - how_many_nt_structure

        if left_or_right == 'left':
            new_structure = structure[:start_position] \
                            + '-' * how_many_nt_structure \
                            + structure[start_position:end_position+1
                                                       +how_many_nt] \
                            + '-' * residual_nt \
                            + structure[end_position+1+how_many_nt:]
        else:
            new_structure = structure[:start_position] \
                            + '-' * residual_nt \
                            + structure[start_position:end_position+1+
                                                       how_many_nt] \
                            + '-' * how_many_nt_structure \
                            + structure[end_position+1+how_many_nt:]
        new_structures.append(new_structure)
    return new_structures


def remove_gaps_same_place(dotbracket_structures):
    if not dotbracket_structures:
        return []
    i = 0
    while i < len(dotbracket_structures[0]):
        all_gap = True
        for structure in dotbracket_structures:
            if structure[i] != '-':
                all_gap = False
                break
        if all_gap:
            new_structures = []
            for structure in dotbracket_structures:
                new_structures.append(structure[:i]+structure[i+1:])
            dotbracket_structures = new_structures
        else:
            i += 1
    return dotbracket_structures


def move_gaps_to_the_loop_centre(dotbracket_structures, structural_blocks):
    new_structures = []
    for dotbracket_structure, blocks in \
            zip(dotbracket_structures, structural_blocks):
        new_structure = list(dotbracket_structure)
        for i in range(len(blocks)-1):
            middle_start_slice = blocks[i][1] + 1
            middle_end_slice = blocks[i+1][0]
            if dotbracket_structure[middle_start_slice -1] == '(' \
                    and dotbracket_structure[middle_end_slice + 1] == ')':
                my_slice = dotbracket_structure[
                           middle_start_slice:middle_end_slice]
                dots = my_slice.count('.')
                gaps = my_slice.count('-')
                new_slice = dots // 2 * '.' + gaps * '-' + \
                            ((dots // 2) + dots % 2) * '.'
                for j in range(len(new_slice)):
                    new_structure[middle_start_slice+j] = new_slice[j]

        new_structures.append(''.join(new_structure))
    # remove - if in all structures
    new_structures = remove_gaps_same_place(new_structures)
    return new_structures


def refine(dotbracket_structures, max_nt, center, repeat):
    for i in range(repeat):
        dotbracket_structures = move_1_2nt_gaps(
            dotbracket_structures, offset=0, max_diff=max_nt, multi_score=1.1)
        dotbracket_structures = remove_gaps_same_place(dotbracket_structures)
        if center:
            representations = [
                structure_to_representation(structure)
                for structure in dotbracket_structures]
            structural_blocks = [
                find_structural_blocks(dotbracket_structure, representation)
                for dotbracket_structure, representation in
                zip(dotbracket_structures, representations)]
            dotbracket_structures = move_gaps_to_the_loop_centre(
                dotbracket_structures, structural_blocks)
    return dotbracket_structures


def refine_from_file(filename, out_filename, max_nt, center, repeat=1):
    file_data = parse_file(filename)
    dotbracket_structures = [x[2] for x in file_data]
    dotbracket_structures = refine(
        dotbracket_structures, max_nt, center, repeat)


    result = convert_to_file_data(file_data, dotbracket_structures)

    with open(out_filename, 'w') as f:
        for element in result:
            f.write("{}\n{}\n{}\n".format(*element))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="Input file (dot bracket)", required=True)
    parser.add_argument("-o", help="Output file (dot bracket)", required=True)
    parser.add_argument(
        "-max_refinement", help="Maximum refinement range (nt)", type=int,
        default=5)
    parser.add_argument(
        '-no_center', help="Center gaps within loops?", action='store_true')
    parser.add_argument(
        '-repeat_refinement', help="How many times to repeat", type=int,
        default=1)
    args = parser.parse_args()
    refine_from_file(
        args.i, args.o, args.max_refinement, not args.no_center,
        args.repeat_refinement)


if __name__ == '__main__':
    main()
