"""FastAPI 应用入口。"""

from fastapi import FastAPI

from .envelope import compute_envelope
from .mixture import candidate_spectrum, fit_proportions, normalize_measured
from .models import (
    CandidateProportion,
    EnvelopeRequest,
    EnvelopeResponse,
    FittedIntensity,
    MixtureFitRequest,
    MixtureFitResponse,
    Peak,
    Residual,
)

app = FastAPI(
    title="Isotope Envelope API",
    summary="同位素峰包络计算与混合谱拟合",
    description=(
        "以各元素最小质量数为基准，按独立抽样卷积生成分子的质量偏移分布；"
        "相同偏移合并概率，并按概率加权给出精确质量质心。"
        "另提供混合谱拟合：在比例非负且和为 1 的约束下，"
        "求使拟合谱与实测谱残差平方和全局最小的候选比例。"
    ),
    version="1.1.0",
)


def element_specs(elements):
    """把元素模型列表转换为 ``(基准质量数加权和, 包络计算参数)``。

    基准 = 各元素最小质量数 × 原子数之和；包络参数中同位素质量数
    已转换为相对各元素最小质量数的偏移。
    """
    base = 0
    specs = []
    for el in elements:
        el_base = min(iso.mass_number for iso in el.isotopes)
        base += el_base * el.count
        isotopes = [
            (iso.mass_number - el_base, iso.exact_mass, iso.abundance)
            for iso in el.isotopes
        ]
        specs.append((el.count, isotopes))
    return base, specs


@app.get("/")
def root():
    return {
        "service": "isotope-envelope",
        "endpoints": ["POST /envelope", "POST /mixture-fit"],
        "docs": "/docs",
    }


@app.post("/envelope", response_model=EnvelopeResponse)
def envelope(request: EnvelopeRequest) -> EnvelopeResponse:
    _, specs = element_specs(request.elements)
    peaks, total = compute_envelope(specs)
    return EnvelopeResponse(
        peaks=[
            Peak(offset=offset, probability=str(prob), centroid=str(centroid))
            for offset, prob, centroid in peaks
        ],
        total_probability=str(total),
    )


@app.post("/mixture-fit", response_model=MixtureFitResponse)
def mixture_fit(request: MixtureFitRequest) -> MixtureFitResponse:
    spectra = []
    for candidate in request.candidates:
        base, specs = element_specs(candidate.elements)
        spectra.append(candidate_spectrum(base, specs))

    measured = normalize_measured(request.mass_number, request.intensity)
    x, grid, fitted, residual, sse = fit_proportions(spectra, measured)

    return MixtureFitResponse(
        proportions=[
            CandidateProportion(id=candidate.id, proportion=float(p))
            for candidate, p in zip(request.candidates, x)
        ],
        fitted=[
            FittedIntensity(mass_number=mass, intensity=float(value))
            for mass, value in zip(grid, fitted)
        ],
        residuals=[
            Residual(mass_number=mass, residual=float(value))
            for mass, value in zip(grid, residual)
        ],
        sse=float(sse),
    )
