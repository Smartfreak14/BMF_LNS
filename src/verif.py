#!/usr/bin/env python3

import argparse
import csv

def read_csv_matrix(filepath):
    matrix = []
    with open(filepath, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            matrix.append([int(float(x)) for x in row if x.strip()])
    return matrix

def boolean_matrix_multiply(A, B):
    
    m = len(A)
    k = len(A[0]) if A else 0
    n = len(B[0]) if B else 0
    
    result = [[0 for _ in range(n)] for _ in range(m)]
    
    for i in range(m):
        for j in range(n):
            for l in range(k):
                if A[i][l] == 1 and B[l][j] == 1:
                    result[i][j] = 1
                    break  
    
    return result

def count_errors(X, AB):
    
    covered_zeros = 0     
    uncovered_ones = 0    
    missing = 0
    m = len(X)
    n = len(X[0]) if X else 0

    for i in range(m):
        for j in range(n):
            if X[i][j] == -1:
                missing += 1
                continue
            if X[i][j] == 0 and AB[i][j] == 1:
                covered_zeros += 1
            elif X[i][j] == 1 and AB[i][j] == 0:
                uncovered_ones += 1

    if missing > 0:
        print(f"(ignored {missing} missing values)")

    total = covered_zeros + uncovered_ones
    return total, covered_zeros, uncovered_ones


parser = argparse.ArgumentParser(description='Count the number of reconstruction errors A 0 B on X')

parser.add_argument('A', metavar='A', type=str, help='path to the matrix A')
parser.add_argument('B', metavar='B', type=str, help='path to the matrix B')
parser.add_argument('X', metavar='X', type=str, help='path to the input matrix X')

args = parser.parse_args()

# Lire les matrices
A = read_csv_matrix(args.A)
B = read_csv_matrix(args.B)
X = read_csv_matrix(args.X)

# Calculer A x B
AB = boolean_matrix_multiply(A, B)

# Compter les erreurs (hard/soft separement)
total, covered_zeros, uncovered_ones = count_errors(X, AB)

k = len(A[0]) if A and A[0] else 0
from_below = "yes" if covered_zeros == 0 else "no"
print(f"k = {k}")
print(f"Covered 0-entries (hard violations) = {covered_zeros}")
print(f"Uncovered 1-entries (soft violations, err_M) = {uncovered_ones}")
print(f"From-below = {from_below}")
print(f"Reconstruction error = {total}")


