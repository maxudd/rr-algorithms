def ore_ord(a):
    try:
        return -1 if a.is_zero() else a.degree()
    except AttributeError:
        return -1 if a == 0 else a.degree()


def row_pivot_info(M, i):
    best_deg = -1
    best_piv = 0
    for j in range(M.ncols()):
        d = ore_ord(M[i, j])
        if d > best_deg:
            best_deg = d
            best_piv = j + 1
    return best_piv, best_deg


def order_matrix(M):
    return matrix(ZZ, M.nrows(), M.ncols(),
                  lambda i, j: ore_ord(M[i, j]))

def pick_qrr_pair_reverse(M):
    buckets = {}
    for i in range(M.nrows()):
        p, deg = row_pivot_info(M, i)
        if p != 0:
            buckets.setdefault(p, []).append((i, deg))

    for p in range(1, M.ncols() + 1):
        if p in buckets and len(buckets[p]) > 1:
            lst = sorted(buckets[p], key=lambda t: (t[1], t[0]))
            (i, deg_i), (j, deg_j) = lst[-1], lst[-2]
            return p, i, deg_i, j, deg_j

    return None


def qrr_one_step(M, U=None):
    choice = pick_qrr_pair_reverse(M)
    if choice is None:
        return False

    p, i, deg_i, j, deg_j = choice
    A = M.base_ring()
    D = A.gen()

    lc_i = M[i, p - 1].leading_coefficient()
    lc_j = M[j, p - 1].leading_coefficient()
    e = deg_i - deg_j

    shift = D**e

    new_row = [
        lc_j * M[i, c] - lc_i * shift * M[j, c]
        for c in range(M.ncols())
    ]
    for c, val in enumerate(new_row):
        M[i, c] = val

    if U is not None:
        new_u_row = [
            lc_j * U[i, c] - lc_i * shift * U[j, c]
            for c in range(U.ncols())
        ]
        for c, val in enumerate(new_u_row):
            U[i, c] = val

    return True

def random_dense_poly(R, d, coeff_bound=5, rng=None):
    if d < 0:
        return R.zero()
    if coeff_bound < 1:
        raise ValueError("coeff_bound must be >= 1")

    if rng is None:
        rng = random.Random()

    x = R.gen()
    coeffs = []
    for _ in range(d + 1):
        c = 0
        while c == 0:
            c = rng.randint(-coeff_bound, coeff_bound)
        coeffs.append(QQ(c))

    return sum(coeffs[k] * x**k for k in range(d + 1))

def pascal_base_rows(R, m, d, beta=None, coeff_bound=5, rng=None):
    if m < 2:
        raise ValueError("m must be >= 2")

    if rng is None:
        rng = random.Random()

    if beta is None:
        beta = m + 1

    common = [R.zero()]
    for j in range(1, m):
        common.append(random_dense_poly(R, d, coeff_bound=coeff_bound, rng=rng))

    rows = []
    for i in range(m):
        row = [R.one()]
        for j in range(1, m):
            row.append(R(binomial(beta + i, j)) + common[j])
        rows.append(row)
    return rows

def reverse_lift_rows(base_rows, A, ell):
    if ell < 0:
        raise ValueError("ell must be >= 0")

    D = A.gen()
    rows = [[A(a) for a in row] for row in base_rows]

    m = len(rows)
    n = len(rows[0])

    for _ in range(ell):
        new_rows = [None] * m

        # b_1^(s+1) = b_1^(s) + D*b_m^(s)
        new_rows[0] = [rows[0][c] + D * rows[m - 1][c] for c in range(n)]

        # b_i^(s+1) = b_i^(s) + b_{i-1}^(s+1), i = 2,...,m
        for i in range(1, m):
            new_rows[i] = [rows[i][c] + new_rows[i - 1][c] for c in range(n)]

        rows = new_rows

    return rows

class Algorithm:
    def __init__(self, name, func):
        self.name = name
        self.func = func

class TestFramework:
    def __init__(self, alg_list):
        self.alg_list = alg_list
        
    def set_test_params(self, N=50):
        self.N = N

    def run_single_test(self, m, ell, d, timeout_seconds):
        while True:
            res = {alg.name: [] for alg in self.alg_list}
            test_mat = self.generate_matrix(A, m, ell, d)
            for alg in self.alg_list:
                start_time = time.time()
                try:
                    alarm(timeout_seconds)
                    _, steps = alg.func(test_mat)
                    cancel_alarm()
                    elapsed = time.time() - start_time
                    res[alg.name].append(elapsed)
                    res[alg.name].append(steps)
                except AlarmInterrupt:
                    cancel_alarm()
                    print('dolgo')
                    break
                except Exception as e:
                    cancel_alarm()
                    print(f"Error: {e}")
                    break
            else:
                return res

    def generate_matrix(self, A, m, ell, d,
                         beta=None,
                         coeff_bound=5,
                         seed=None,
                         max_tries=20):

        if m < 2:
            raise ValueError("m must be >= 2")
    
        rng = random.Random(seed)
    
        for _ in range(max_tries):
            base_rows = pascal_base_rows(
                A.base_ring(), m, d,
                beta=beta,
                coeff_bound=coeff_bound,
                rng=rng
            )
    
            lifted_rows = reverse_lift_rows(base_rows, A, ell)
            L = matrix(A, lifted_rows)
    
            return L
    
    def run_experiment(self, N, param_name, param_values, fixed_params, timeout=60):
        results = {
            alg.name: {
                "time": [],
                "steps": []
            } for alg in self.alg_list
        }
    
        for val in param_values:
            print(f"=== {param_name}: {val} ===")
    
            agg = {
                alg.name: {"time": [], "steps": []}
                for alg in self.alg_list
            }
    
            for t in range(N):
                params = fixed_params.copy()
                params[param_name] = val
    
                m = params.get("m")
                ell = params.get("ell")
                d = params.get("d")
                
                single_res = self.run_single_test(m, ell, d, timeout)
                for alg in self.alg_list:
                    agg[alg.name]["time"].append(single_res[alg.name][0])
                    agg[alg.name]["steps"].append(single_res[alg.name][1])
    
            for alg in self.alg_list:
                times = agg[alg.name]["time"]
                steps = agg[alg.name]["steps"]
    
                results[alg.name]["time"].append({
                    "avg": np.mean(times),
                    "max": np.max(times),
                    "min": np.min(times)
                })
    
                results[alg.name]["steps"].append({
                    "avg": np.mean(steps),
                    "max": np.max(steps)
                })
    
        return results
    
    def plot_results(self, param_name, param_values, results, log=False):
        fig, ax = plt.subplots()
    
        if log:
            ax.set_yscale("log")
    
        for alg_name, data in results.items():
            avg = [x["avg"] for x in data["time"]]
            # mx  = [x["max"] for x in data["time"]]
            # mn  = [x["min"] for x in data["time"]]
    
            ax.plot(param_values, avg, ".-", label=f"{alg_name} avg")
            # ax.plot(param_values, mx, "--", label=f"{alg_name} max")
            # ax.plot(param_values, mn, ":", label=f"{alg_name} min")
    
        ax.set_xlabel(param_name)
        ax.set_ylabel("Время (с)")
        ax.set_xticks(param_values)
        ax.grid()
        ax.legend()
        plt.show()

    def plot_results_separate(self, param_name, param_values, results, log=False):
        for alg_name, data in results.items():
            fig, ax = plt.subplots()
    
            if log:
                ax.set_yscale("log")
    
            avg = [x["avg"] for x in data["time"]]
            mx  = [x["max"] for x in data["time"]]
            mn  = [x["min"] for x in data["time"]]
    
            ax.plot(param_values, avg, ".-", color="b", label=f"{alg_name} avg")
            ax.plot(param_values, mx, ".-", color="r", label=f"{alg_name} max")
            ax.plot(param_values, mn, ".-", color="g", label=f"{alg_name} min")
    
            ax.set_xlabel(param_name)
            ax.set_ylabel("Время (с)")
            ax.set_xticks(param_values)
            ax.grid()
            ax.legend()
            plt.show()

    def save_results(self, results, filename):
        with open(filename, "wb") as f:
            pickle.dump(results, f)

    def load_results(self, filename):
        with open(filename, "rb") as f:
            return pickle.load(f)