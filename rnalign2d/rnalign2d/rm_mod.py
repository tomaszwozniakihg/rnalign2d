import argparse
try:
    from .rnalign2d import remove_modifications
except SystemError:
    from rnalign2d.rnalign2d import remove_modifications


def unmodify_file(filename, out_filename):
    result = []
    name = None
    sequence = None
    structure = None
    counter = 0
    with open(filename, 'r') as f:
        for line in f.readlines():
            if counter % 3 == 0:
                if counter != 0 and len(line.strip()) > 0:
                    result.append(name)
                    result.append(remove_modifications(sequence))
                    result.append(structure)
                name = line.strip()
            elif counter % 3 == 1:
                sequence = line.strip()
            else:
                structure = line.strip()
            counter += 1
        result.append(name)
        result.append(remove_modifications(sequence))
        result.append(structure)
    with open(out_filename, 'w') as f:
        f.write('\n'.join(result))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="Input file (dot bracket)",
                        required=True)
    parser.add_argument("-o", help="Output file (dot bracket)",
                        required=True)
    args = parser.parse_args()
    unmodify_file(args.i, args.o)


if __name__ == '__main__':
    main()
