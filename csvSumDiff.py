#! /usr/bin/python3

import sys
import os
from math import sqrt

def values(file_name):
    with open(file_name, "r") as f:
        for l in f.readlines():
            for v in (s for s in l.split(",") if s.strip()!=""):
                try:
                    yield float(v)
                except ValueError:
                    print("could not convert '"+str(v)+"' to float")

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
    variance = total_sq/n - total*total/n/n
    print("total=", total)
    print("avg=", total/n)
    print("variance=", variance)
    print("std-dev=", sqrt(variance))


if __name__ == "__main__":
    main()
