#! /usr/bin/python3

import sys
import os

def values(file_name):
    with open(file_name, "r") as f:
        for l in f.readlines():
            for v in l.split(","):
                try:
                    yield float(v)
                except:
                    pass

def main():
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    print("comparing "+file1+" and "+file2)
    total = 0
    total_sq = 0
    n = 0
    for v1, v2 in zip(values(file1), values(file2)):
        diff = v2-v1
        if diff*diff > 1:
            print(v1,v2)
        total += diff
        total_sq += diff*diff
        n += 1
    print("total=", total)
    print("avg=", total/n)
    print("variance=", n*total_sq/total/total)


if __name__ == "__main__":
    main()
