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

def pad_None(iterator):
    for v in iterator:
        yield v
    yield None

def main():
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    print("comparing "+file1+" and "+file2+"\n")
    total = 0
    total_sq = 0
    n = 0
    for v1, v2 in zip(pad_None(values(file1)), pad_None(values(file2))):
        if v1==None or v2==None:
            if v1!=None or v2!=None:
                print("\033[1;31mFiles aren't the same length\033[0m\n")
            break
        diff = v2-v1
        if diff*diff > 1:
            print("\033[31m"+(str(v1)+",").ljust(8)+"\t"+str(v2)+"\033[0m")
        elif diff*diff > 0.01:
            print("\033[33m"+(str(v1)+",").ljust(8)+"\t"+str(v2)+"\033[0m")
        elif diff*diff > 0.0000001:
            print((str(v1)+",").ljust(8)+"\t"+str(v2))
        total += diff
        total_sq += diff*diff
        n += 1
    variance = total_sq/n - total*total/n/n
    print("\ntotal=", total)
    print("avg=", total/n)
    print("variance=", variance)
    print("std-dev=", sqrt(variance))


if __name__ == "__main__":
    main()
