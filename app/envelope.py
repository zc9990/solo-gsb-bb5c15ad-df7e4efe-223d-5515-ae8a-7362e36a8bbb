"""同位素峰包络核心计算。

以各元素最小质量数为基准（质量偏移 0），把每个原子看作独立抽样，
通过卷积得到整个分子的质量偏移分布。

分布用 ``dict[offset, (probability, mass_moment)]`` 表示：

- ``probability``：落在该偏移上的总概率（相同偏移合并）。
- ``mass_moment``：``probability * 这些同位素组合的精确质量之和``，
  用于在最后按概率加权求精确质量质心。

全程使用高精度 ``Decimal``，不裁剪任何低概率中间项。
"""

from decimal import Decimal, localcontext

# 计算精度（有效数字位数）。60 个原子、丰度最极端时概率可达 1e-600 量级，
# Decimal 的指数范围足够大，80 位有效数字可保证质心与概率和的精度。
PRECISION = 80

# 分布类型：offset -> (probability, mass_moment)
Distribution = dict[int, tuple[Decimal, Decimal]]


def convolve(d1: Distribution, d2: Distribution) -> Distribution:
    """卷积两个分布（对应两组相互独立的原子）。"""
    out: Distribution = {}
    for off1, (p1, m1) in d1.items():
        if p1 == 0:
            continue
        for off2, (p2, m2) in d2.items():
            if p2 == 0:
                continue
            key = off1 + off2
            prob = p1 * p2
            moment = m1 * p2 + p1 * m2
            if key in out:
                acc_p, acc_m = out[key]
                out[key] = (acc_p + prob, acc_m + moment)
            else:
                out[key] = (prob, moment)
    return out


def element_distribution(count: int, isotopes: list[tuple[int, Decimal, Decimal]]) -> Distribution:
    """单个元素的分布。

    ``isotopes`` 为 ``(质量偏移, 精确质量, 丰度)`` 列表，
    质量偏移已以该元素最小质量数为基准。``count`` 个原子逐一卷积（独立抽样）。
    """
    dist: Distribution = {0: (Decimal(1), Decimal(0))}
    single: Distribution = {}
    for offset, exact_mass, abundance in isotopes:
        if abundance == 0:
            continue  # 概率严格为 0 的项对结果无任何贡献
        single[offset] = (abundance, abundance * exact_mass)
    for _ in range(count):
        dist = convolve(dist, single)
    return dist


def compute_envelope(elements: list[tuple[int, list[tuple[int, Decimal, Decimal]]]]):
    """计算分子整体的峰包络。

    参数：``[(原子数, [(质量偏移, 精确质量, 丰度), ...]), ...]``。
    返回：``(peaks, total_probability)``，其中 ``peaks`` 为按偏移升序的
    ``(offset, probability, centroid)`` 列表，仅含非零概率。
    """
    with localcontext() as ctx:
        ctx.prec = PRECISION

        molecule: Distribution = {0: (Decimal(1), Decimal(0))}
        for count, isotopes in elements:
            if count == 0:
                continue
            molecule = convolve(molecule, element_distribution(count, isotopes))

        peaks = []
        total = Decimal(0)
        for offset in sorted(molecule):
            prob, moment = molecule[offset]
            if prob == 0:
                continue
            peaks.append((offset, prob, moment / prob))
            total += prob
        return peaks, total
