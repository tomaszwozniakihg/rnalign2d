#!/usr/bin/env python
import argparse
from collections import defaultdict
try:
    from .conversion import SIMPLE_CONVERSION, PSEUDOKNOT_CONVERSION, LETTERS
except SystemError:
    from rnalign2d.conversion import SIMPLE_CONVERSION, PSEUDOKNOT_CONVERSION, \
        LETTERS


OPENING_BRACKETS = "([{<ABCD"
CLOSING_BRACKETS = ")]}>abcd"


def score_brackets(
        dot_bracket1, dot_bracket2, score_same_brackets, score_other_brackets,
        score_reverse_brackets, score_brackets_dots, score_two_dots):
    score = 0
    if dot_bracket1 == dot_bracket2 and dot_bracket1 == '.':
        score = score_two_dots
    if dot_bracket1 == dot_bracket2:
        score = score_same_brackets
    elif (dot_bracket1 in OPENING_BRACKETS and
                  dot_bracket2 in OPENING_BRACKETS) or (
                    dot_bracket1 in CLOSING_BRACKETS and
                    dot_bracket2 in CLOSING_BRACKETS):
        score = score_other_brackets
    elif (dot_bracket1 in OPENING_BRACKETS and
                  dot_bracket2 in CLOSING_BRACKETS) or (
                    dot_bracket1 in CLOSING_BRACKETS and
                    dot_bracket2 in OPENING_BRACKETS):
        score = score_reverse_brackets
    elif dot_bracket1 == '.' or dot_bracket2 == '.':
        score = score_brackets_dots
    return score


def create_matrix(
        score_same_brackets, score_other_brackets, score_reverse_brackets,
        score_brackets_dots, score_two_dots, add_score_for_seq_match,
        mode='simple'):
    """
    Function that create matrix that can be used for further analysis,
    please take note, that mode must be the same in case of matrix and
    multiple sequence alignment, otherwise random-like effects will occur

    :param score_same_brackets: int, score for the same tye of brackets
    like ( and (
    :param score_other_brackets: int, score for different type of brackets
    like ( and [
    :param score_reverse_brackets: int, socre for reverse brackets like ( and )
    :param score_brackets_dots: int, socre for brakcet and dot like ( and .
    :param score_two_dots: int, score for two dots like . and .
    :param add_score_for_seq_match: int, value to add
    if sequence letter is the same
    :param mode: string, simple - only level one pseudoknots, pseudo -
    multiplelevel of pseudoknots
    :return: string containing matrix that can be saved
    """
    header = "       A    C    D    E    F    G    H    I    K    L    M    " \
             "N    P    Q    R    S    T    V    W    Y"
    matrix = defaultdict(dict)
    if mode == 'simple':
        for letter1 in LETTERS:
            nucleotide1 = None
            dot_bracket1 = None
            for nucleotide in SIMPLE_CONVERSION:
                for dot_bracket in SIMPLE_CONVERSION[nucleotide]:
                    if SIMPLE_CONVERSION[nucleotide][dot_bracket] == letter1:
                        nucleotide1 = nucleotide
                        dot_bracket1 = dot_bracket

            for letter2 in LETTERS:
                nucleotide2 = None
                dot_bracket2 = None
                for nucleotide in SIMPLE_CONVERSION:
                    for dot_bracket in SIMPLE_CONVERSION[nucleotide]:
                        if SIMPLE_CONVERSION[nucleotide][dot_bracket] == \
                                letter2:
                            nucleotide2 = nucleotide
                            dot_bracket2 = dot_bracket
                score = score_brackets(
                    dot_bracket1, dot_bracket2, score_same_brackets,
                    score_other_brackets, score_reverse_brackets,
                    score_brackets_dots, score_two_dots)
                if nucleotide1 == nucleotide2:
                    score += add_score_for_seq_match
                matrix[letter1][letter2] = score
    elif mode == 'pseudo':
        for letter1 in LETTERS:
            dot_bracket1 = None
            for dot_bracket in PSEUDOKNOT_CONVERSION:
                if PSEUDOKNOT_CONVERSION[dot_bracket] == letter1:
                    dot_bracket1 = dot_bracket
            for letter2 in LETTERS:
                score = 0
                dot_bracket2 = None
                for dot_bracket in PSEUDOKNOT_CONVERSION:
                    if PSEUDOKNOT_CONVERSION[dot_bracket] == letter2:
                        dot_bracket2 = dot_bracket
                if dot_bracket2 is not None and dot_bracket1 is not None:
                    score = score_brackets(
                        dot_bracket1, dot_bracket2, score_same_brackets,
                        score_other_brackets, score_reverse_brackets,
                        score_brackets_dots, score_two_dots)
                matrix[letter1][letter2] = score
    else:
        print('Wrong mode')
    text = [header]
    for letter1 in LETTERS:
        string = [letter1, '  ']
        for letter2 in LETTERS:
            score = matrix[letter1][letter2]
            string.append(str(score).rjust(5))
        text.append("".join(string))
    return "\n".join(text)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-same", type=int, help="Score same brackets like ( and (", default=7)
    parser.add_argument(
        "-other", type=int, help="Score other brackets like [ and (", default=2)
    parser.add_argument(
        "-reverse", type=int, help="Score reverse brackets like ) and (",
        default=-10)
    parser.add_argument(
        "-bracket_dot", type=int, help="Score bracket and dot like . and (",
        default=-1)
    parser.add_argument(
        "-dot_dot", type=int, help="Score two dota like . and .", default=3)
    parser.add_argument(
        "-seq_match_add", type=int, help="Add to score for sequence match",
        default=1)
    parser.add_argument(
        "-mode", help="Mode for matrix creation - it can be 'simple' for "
                      "comparision with maximum of one level of pseudoknots "
                      "or 'pseudo' for multiple level of pseudoknots",
        choices=['simple', 'pseudo'], default='simple')
    parser.add_argument("-o", help="Output file", default='result_matrix')
    args = parser.parse_args()
    data = create_matrix(
        args.same, args.other, args.reverse, args.bracket_dot, args.dot_dot,
        args.seq_match_add, args.mode)
    with open(args.o, 'w') as f:
        f.write(data)


if __name__ == '__main__':
    main()
