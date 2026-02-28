#!/usr/bin/env python3
"""
Vérification de la factorisation BMF: A ◦ B ≈ X
Calcule le nombre d'erreurs de reconstruction.

Version sans dépendances externes (numpy/pandas).
"""

import argparse
import csv

def read_csv_matrix(filepath):
    """Lit une matrice depuis un fichier CSV."""
    matrix = []
    with open(filepath, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            matrix.append([int(float(x)) for x in row if x.strip()])
    return matrix

def boolean_matrix_multiply(A, B):
    """
    Calcule le produit booléen A ◦ B.
    (A ◦ B)[i,j] = OR_{l} (A[i,l] AND B[l,j])
    """
    m = len(A)
    k = len(A[0]) if A else 0
    n = len(B[0]) if B else 0
    
    result = [[0 for _ in range(n)] for _ in range(m)]
    
    for i in range(m):
        for j in range(n):
            for l in range(k):
                if A[i][l] == 1 and B[l][j] == 1:
                    result[i][j] = 1
                    break  # OR: dès qu'on trouve un 1, on peut arrêter
    
    return result

def count_errors(X, AB):
    """
    Compte le nombre de cellules différentes entre X et A◦B.
    Ignore les valeurs manquantes (X[i][j] == -1).
    """
    errors = 0
    missing = 0
    m = len(X)
    n = len(X[0]) if X else 0
    
    for i in range(m):
        for j in range(n):
            if X[i][j] == -1:  # Valeur manquante → ignorer
                missing += 1
                continue
            if X[i][j] != AB[i][j]:
                errors += 1
    
    if missing > 0:
        print(f"(ignored {missing} missing values)")
    
    return errors


parser = argparse.ArgumentParser(description='Count the number of reconstruction errors A ◦ B on X')

parser.add_argument('A', metavar='A', type=str, help='path to the matrix A')
parser.add_argument('B', metavar='B', type=str, help='path to the matrix B')
parser.add_argument('X', metavar='X', type=str, help='path to the input matrix X')

args = parser.parse_args()

# Lire les matrices
A = read_csv_matrix(args.A)
B = read_csv_matrix(args.B)
X = read_csv_matrix(args.X)

# Calculer A ◦ B
AB = boolean_matrix_multiply(A, B)

# Compter les erreurs
error = count_errors(X, AB)

k = len(A[0]) if A and A[0] else 0
print(f"k = {k}")
print(f"Reconstruction error = {error}")


