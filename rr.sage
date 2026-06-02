def null_space_vector(M, nonzero=False):
    if nonzero:
        nonzero_indices = [i for i in range(M.nrows()) if not M.row(i).is_zero()]
        M_nonzero = M.matrix_from_rows(nonzero_indices)
        k = kernel(M_nonzero)
    else:
        k = kernel(M)
    return list(k.basis_matrix()[0]) if k else []

def find_t_index(v, delta, n):
    appr_delta = [-1 if v[i] == 0 else delta[i] for i in range(n)]

    return appr_delta.index(max(appr_delta))

def frMat_and_rowRanks(M):
    rrs = []
    fm = []
    for row in M:
        rr = max(map(lambda e: 0 if e == 0 else e.degree(), row))
        rrs.append(rr)
        fr = list(map(lambda e: e.leading_coefficient() if e.degree() == rr else 0, row))
        fm.append(fr)
    return (matrix(fm), rrs)

def row_red(M, transformation=False):
    n = M.nrows()
    L = copy(M)
    # for i in range(n):
    #    reduce_row_coefficients(L, i)
    if transformation: 
        U = identity_matrix(A, n)
    frM, rr = frMat_and_rowRanks(L)
    step = 0
    while null_space_vector(frM, nonzero=True):
        # print(f'=======Step {step}=========')
        step += 1
        v = null_space_vector(frM)
        t = find_t_index(v, rr, n)
        # print(rr)
        # print(f'v: {v}')
        L_t_th_row = sum ([ (0 if v[i] == 0 else v[i]*D^(rr[t]-rr[i])) * L[i] for i in range(n)])
        if transformation:
            U_t_th_row = sum ([ (0 if v[i] == 0 else v[i]*D^(rr[t]-rr[i])) * U[i] for i in range(n)])
            U.set_row(t,U_t_th_row)
        L.set_row(t,L_t_th_row)
        # reduce_row_coefficients(L, t)
        frM, rr = frMat_and_rowRanks(L)
    return (L, U, step) if transformation else (L, step)
