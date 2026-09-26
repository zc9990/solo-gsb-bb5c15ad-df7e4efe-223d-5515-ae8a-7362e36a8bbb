"""混合谱拟合核心计算。

在候选谱与实测谱的质量数并集上补零后求解：

    min  ||A x - b||^2
    s.t. x >= 0 且 sum(x) == 1

其中 A 的各列为候选的归一化概率谱，b 为实测归一化谱；多解时取
比例平方和 ||x||^2 最小者。

候选数仅 2~8，枚举全部非空子集（至多 255 个）：对每个子集解等式约束
（比例和为 1）最小二乘的 KKT 系统——对称、可能奇异——用 Jacobi 特征
分解取最小范数解，并检验非负性。由于最优拟合向量 y* = A x 唯一，全局
最优解集是凸集 {x : A x = y*, sum(x) = 1, x >= 0}，其最小范数点必在其
支撑集对应的仿射子空间上取得，因此枚举法能精确找到"先最小 RSS、再
最小 ||x||"的字典序最优解。

全程 float64；谱与比例均为 O(1) 量级，数值误差远低于 1e-8。
无第三方依赖，仅使用标准库。
"""

from __future__ import annotations

from itertools import combinations

# 非负可行性容差：float64 噪声约 1e-16，比例真值为 0 时解出值在 0 附近。
FEAS_TOL = 1e-9
# RSS 平局判定（相对）：远大于格拉姆矩阵点积的累计舍入（~1e-14）。
RSS_TOL = 1e-12
# 最终比例裁剪阈值：|x| 小于它视为 0。
CLIP_TOL = 1e-12
# Jacobi 特征分解的零特征值判定（相对最大特征值）。
EIG_ZERO_TOL = 1e-12


def _jacobi_eigh(matrix: list[list[float]]) -> tuple[list[float], list[list[float]]]:
    """对称矩阵特征分解（经典 Jacobi 旋转，每次消去最大非对角元）。

    返回 ``(特征值列表, 特征向量按列存储的矩阵)``。矩阵规模不超过 9，
    二次收敛，特征值精度约 1e-15 相对量级。
    """
    n = len(matrix)
    a = [row[:] for row in matrix]
    v = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for _ in range(1000):
        p, q, off = 0, 0, 0.0
        for i in range(n):
            for j in range(i + 1, n):
                if abs(a[i][j]) > off:
                    p, q, off = i, j, abs(a[i][j])
        if off <= 1e-15:
            break
        theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q])
        t = (1.0 if theta >= 0 else -1.0) / (abs(theta) + (theta * theta + 1.0) ** 0.5)
        c = 1.0 / (t * t + 1.0) ** 0.5
        s = t * c
        for k in range(n):  # a = a @ J（列变换）
            akp, akq = a[k][p], a[k][q]
            a[k][p] = c * akp - s * akq
            a[k][q] = s * akp + c * akq
        for k in range(n):  # a = J^T @ a（行变换）
            apk, aqk = a[p][k], a[q][k]
            a[p][k] = c * apk - s * aqk
            a[q][k] = s * apk + c * aqk
        for k in range(n):  # v = v @ J（累计特征向量）
            vkp, vkq = v[k][p], v[k][q]
            v[k][p] = c * vkp - s * vkq
            v[k][q] = s * vkp + c * vkq
    return [a[i][i] for i in range(n)], v


def _min_norm_solve(kkt: list[list[float]], rhs: list[float]) -> list[float]:
    """对称（可奇异）系统 ``kkt @ z = rhs`` 的最小范数最小二乘解。

    对特征值按 ``EIG_ZERO_TOL`` 截断（ Moore-Penrose 伪逆），
    奇异时返回所有解中范数最小者。
    """
    eigvals, vecs = _jacobi_eigh(kkt)
    lam_max = max(abs(lam) for lam in eigvals)
    zero_tol = EIG_ZERO_TOL * lam_max
    n = len(rhs)
    z = [0.0] * n
    for i in range(n):
        lam = eigvals[i]
        if abs(lam) <= zero_tol:
            continue
        coeff = sum(vecs[r][i] * rhs[r] for r in range(n)) / lam
        for r in range(n):
            z[r] += coeff * vecs[r][i]
    return z


def fit_proportions(columns: list[list[float]], b: list[float]) -> list[float]:
    """求非负且和为 1 的比例，使拟合谱与实测谱的残差平方和全局最小。

    ``columns``：各候选在并集网格上的归一化谱（k 个长度 m 的列表）；
    ``b``：实测归一化谱（长度 m）。多解时返回比例平方和最小者。
    返回长度 k 的比例列表（和为 1，非负，|x| < CLIP_TOL 的项裁为 0）。
    """
    k = len(columns)
    m = len(b)

    # 一次性计算完整格拉姆矩阵 G = A^T A、 c = A^T b、bb = b·b；
    # 各子集问题只需取对应子矩阵 / 子向量，RSS 用二次型求值。
    g = [
        [sum(columns[i][r] * columns[j][r] for r in range(m)) for j in range(k)]
        for i in range(k)
    ]
    c = [sum(columns[i][r] * b[r] for r in range(m)) for i in range(k)]
    bb = sum(v * v for v in b)

    feasible = []  # (rss, norm2, subset, x)
    for size in range(1, k + 1):
        for subset in combinations(range(k), size):
            # 等式约束最小二乘的 KKT 系统 [G_S 1; 1^T 0] [x; λ] = [c_S; 1]
            n = size + 1
            kkt = [[0.0] * n for _ in range(n)]
            for ai, i in enumerate(subset):
                for aj, j in enumerate(subset):
                    kkt[ai][aj] = g[i][j]
                kkt[ai][size] = 1.0
                kkt[size][ai] = 1.0
            rhs = [c[i] for i in subset] + [1.0]
            x = _min_norm_solve(kkt, rhs)[:size]
            if any(v < -FEAS_TOL for v in x):
                continue
            # rss = x^T G_S x - 2 c_S·x + bb（浮点误差下钳到非负）
            rss = bb
            for ai, i in enumerate(subset):
                rss -= 2.0 * c[i] * x[ai]
                for aj, j in enumerate(subset):
                    rss += x[ai] * g[i][j] * x[aj]
            if rss < 0.0:
                rss = 0.0
            norm2 = sum(v * v for v in x)
            feasible.append((rss, norm2, subset, x))

    # 先取 RSS 全局最小，再在平局（容差内）中取 ||x||^2 最小者
    r_star = min(f[0] for f in feasible)
    tol = RSS_TOL * max(1.0, r_star)
    _, _, subset, x = min(
        (f for f in feasible if f[0] <= r_star + tol), key=lambda f: f[1]
    )

    result = [0.0] * k
    for idx, v in zip(subset, x):
        result[idx] = 0.0 if abs(v) < CLIP_TOL else v
    return result
