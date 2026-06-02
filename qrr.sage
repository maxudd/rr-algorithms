from collections import defaultdict 
from sage.arith.misc import gcd

def reduce_row_coefficients(M, i):
    row_coeffs = []
    for p in M.row(i):
        # p.coefficients() returns a list of coeffs (sparse or dense, depending on your ring)
        row_coeffs.extend([c for c in p.coefficients() if c != 0])  # Use abs to handle negatives; skip zeros
    # print(row_coeffs)
    g = gcd(row_coeffs) if row_coeffs else 0
    # print(g)
    
    if g > 1:
        # Now divide each polynomial's coefficients by g
        for j in range(M.ncols()):
            # print(f'col {j}')
            p = M[i, j]
            if p != 0:
                # map_coefficients applies a function to each coeff
                coeffs = p.list()
                # print(coeffs)
                new_coeffs = [c // g for c in coeffs]
                M[i, j] = A(new_coeffs)

true_degree = lambda elem: elem.degree() if elem else -1

def find_pivots(M, n, num, ps):
    pivot_i = -1
    pivot_d = -1
    for i in range(n):
        d = true_degree(M[num][i])
        if d > pivot_d:
            pivot_i = i
            pivot_d = d
    ps[pivot_i].append((num, pivot_d))

def quick_row_red(N, method='division-free', transformation=False):
    M = copy(N)
    m = M.nrows()
    n = M.ncols()
    for i in range(m):
        reduce_row_coefficients(M, i)
    if transformation:
        U = identity_matrix(A, m)
    pivots = defaultdict(list)
    for i in range(m):
        find_pivots(M, n, i, pivots)
    piv_i = 0
    steps = 0
    while 1:
        if len(pivots[piv_i]) > 1:
            # print(f'=======Step {steps}=========')
            steps += 1
            i,deg_i = pivots[piv_i].pop()
            j,deg_j = pivots[piv_i].pop()
            if deg_j > deg_i:
                i,j = j,i
                deg_i,deg_j = deg_j,deg_i
            # print(i, j)
            c1 = M[i][piv_i].leading_coefficient()
            c2 = M[j][piv_i].leading_coefficient()
            e = deg_i - deg_j
            if method == 'division-free':
                M.set_row(i, c2*M[i])
                M.add_multiple_of_row(i, j, -c1*D^e, 0)
                reduce_row_coefficients(M, i)
                if transformation:
                    U.set_row(i, c2*U[i])
                    U.add_multiple_of_row(i, j, -c1*D^e, 0) 
            else:
                c = - c1 / c2
                coeff = c * D^e
                M.add_multiple_of_row(i, j, coeff, 0)
                if transformation:
                    U.add_multiple_of_row(i, j, coeff, 0)    
            find_pivots(M, n, i, pivots)
            pivots[piv_i].append((j,deg_j))
            piv_i = 0
        else:
            if piv_i == n:
                return M, steps if not transformation else (M, U, steps)
            else:
                piv_i += 1