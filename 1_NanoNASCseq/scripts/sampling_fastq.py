#!/usr/bin/env python
import sys
import random

def main():
    p, r = sys.argv[1:]
    p, r = float(p), int(r)
    
    random.seed(p + r)
    
    lines = list()
    for i, line in enumerate(sys.stdin):
        j = i % 4
        lines.append(line)
        if j == 3:
            if random.random() < p:
                for x in lines:
                    sys.stdout.write(x)
            lines = list()
            

if __name__ == "__main__":
    main()