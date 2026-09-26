"""同位素混合谱拟合核心计算。

给定 2~8 个候选分子（各自的同位素组成）与一条实测谱（整数质量数 -> 强度），
求非负且和为 1 的候选比例，使拟合谱与实测谱的残差平方和（SSE）全局最小。

流程：

1. 复用 :func:`app.envelope.compute_envelope` 得到每个候选的质量偏移分布，
   以「各元素最小质量数 × 原子数之和 + 偏移」转换为绝对整数质量数谱，
   并按概率总和归一化；实测强度同样归一化，双方在质量数并集上补零对齐；
2. 解凸二次规划 ``min ||A x - y||^2``，约束 ``x >= 0``、``sum(x) = 1``。
   候选数 n <= 8，枚举全部 2^n - 1 个非空支撑集，分别解对应的等式约束
   KKT 系统并检验完整 KKT 条件——凸问题中满足 KKT 条件即全局最优；
   最优解不唯一（退化）时，取比例平方和最小的解。
"""

from decimal import Decimal, localcontext
from itertools import combinations

import numpy as np

from .envelope import PRECISION, compute_envelope

# 完整 KKT 检验容差（平稳性 / 对偶可行性）
_KKT_TOL = 1e-7
# KKT 线性系统残差容差
_SOLVE_TOL = 1e-9
# 比例非负可行性容差
_FEAS_TOL = 1e-10
# 目标值并列（视为多解）的容差：差值不超过它时取比例平方和最小者
_TIE_TOL = 1e-12
# 退化情形下重试用的微小正则化系数
_REG = 1e-10


def candidate_spectrum(base: int, elements: list) -> dict[int, float]:
    """单个候选的归一化谱：绝对整数质量数 -> 概率。

    ``base`` 为各元素最小质量数乘原子数之和；``elements`` 为
    ``compute_envelope`` 所需的 ``[(原子数, [(偏移, 精确质量, 丰度)])]``。
    概率按总和归一化（理论上总和即为 1）。
    """
    peaks, total = compute_envelope(elements)
    with localcontext() as ctx:
        ctx.prec = PRECISION
        return {base + offset: float(prob / total) for offset, prob, _ in peaks}


def normalize_measured(mass_numbers: list[int], intensities: list[Decimal]) -> dict[int, float]:
    """实测谱归一化：整数质量数 -> 占总强度的比例。"""
    with localcontext() as ctx:
        ctx.prec = PRECISION
        total = sum(intensities, Decimal(0))
        return {m: float(i / total) for m, i in zip(mass_numbers, intensities)}


def _solve_subset(H: np.ndarray, b: np.ndarray, subset: tuple[int, ...], reg: float):
    """解指定支撑集的等式约束 KKT 系统。

    返回通过完整 KKT 检验的完整比例向量，否则返回 ``None``。
    """
    n = len(b)
    k = len(subset)
    idx = list(subset)
    G = H[np.ix_(idx, idx)]
    if reg:
        G = G + reg * np.eye(k)
    # KKT 系统：[G 1; 1^T 0] [x_S; λ] = [b_S; 1]
    K = np.zeros((k + 1, k + 1))
    K[:k, :k] = G
    K[:k, k] = 1.0
    K[k, :k] = 1.0
    rhs = np.empty(k + 1)
    rhs[:k] = b[idx]
    rhs[k] = 1.0
    # lstsq：奇异（退化）时取最小范数解，正好对应多解时的取舍规则
    sol = np.linalg.lstsq(K, rhs, rcond=None)[0]
    if np.max(np.abs(K @ sol - rhs)) > _SOLVE_TOL:
        return None
    x_s = sol[:k]
    if x_s.min() < -_FEAS_TOL:
        return None
    x = np.zeros(n)
    x[idx] = np.clip(x_s, 0.0, None)
    total = x.sum()
    if total <= 0:
        return None
    x /= total
    # 完整 KKT 条件检验（凸问题 <=> 全局最优）
    g = H @ x - b
    support = x > 0
    mu = np.mean(g[support])
    if np.max(np.abs(g[support] - mu)) > _KKT_TOL:
        return None
    if np.any(g[~support] < mu - _KKT_TOL):
        return None
    return x


def fit_proportions(spectra: list[dict[int, float]], measured: dict[int, float]):
    """求最优候选比例。

    ``spectra``：各候选的 ``{质量数: 概率}``；``measured``：``{质量数: 强度比例}``。
    返回 ``(x, grid, fitted, residual, sse)``：候选比例、升序质量数网格、
    拟合强度、实测减拟合的残差、残差平方和。
    """
    n = len(spectra)
    grid = sorted(set(measured).union(*(s.keys() for s in spectra)))
    A = np.zeros((len(grid), n))
    for j, spec in enumerate(spectra):
        for i, mass in enumerate(grid):
            A[i, j] = spec.get(mass, 0.0)
    y = np.array([measured.get(mass, 0.0) for mass in grid])

    H = A.T @ A
    b = A.T @ y

    # 枚举全部非空支撑集：先精确求解，仅当全部失败（极端数值奇异）时
    # 再以微小正则化重试。汇集所有满足 KKT 条件的候选点——凸问题中它们都是
    # 全局最优。
    pool = []
    for k in range(1, n + 1):
        for subset in combinations(range(n), k):
            x = _solve_subset(H, b, subset, 0.0)
            if x is not None:
                r = A @ x - y
                pool.append((float(r @ r), float(x @ x), x))
    if not pool:
        for k in range(1, n + 1):
            for subset in combinations(range(n), k):
                x = _solve_subset(H, b, subset, _REG)
                if x is not None:
                    r = A @ x - y
                    pool.append((float(r @ r), float(x @ x), x))
    if not pool:  # 理论上不可达：凸二次规划必存在可解的 KKT 点
        raise RuntimeError("混合拟合未能找到可行解")

    sse_min = min(p[0] for p in pool)
    tied = [p for p in pool if p[0] <= sse_min + _TIE_TOL]
    # 多解时取比例平方和最小者
    _, _, x = min(tied, key=lambda p: p[1])

    # 数值噪声（约机器精度）产生的极小非零分量清零，避免退化情形下
    # 这种伪分量因平方和略小而被误选为「最小范数解」；量级远小于 1e-8。
    x = np.where(x < 1e-12, 0.0, x)
    x /= x.sum()
    fitted = A @ x
    residual = y - fitted
    sse = float(residual @ residual)
    return x, grid, fitted, residual, sse
