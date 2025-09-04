#!/usr/bin/env python
import sys

MAPPER = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def main():
    name, seq, qua = None, None, None

    for i, line in enumerate(sys.stdin):
        j = i % 4
        if j == 0:
            name = line[:-1]
        elif j == 1:
            seq = line[:-1]
        elif j == 3:
            qua = line[:-1]

            seq = "".join([MAPPER[b] for b in seq[::-1]])
            qua = qua[::-1]
            sys.stdout.write("%s\n%s\n+\n%s\n" % (name, seq, qua))

        
if __name__ == '__main__':
    main()
    