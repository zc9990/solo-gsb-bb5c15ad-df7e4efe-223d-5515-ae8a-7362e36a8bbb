"""FastAPI 应用入口。"""

from fastapi import FastAPI

from .envelope import compute_envelope
from .mixture import fit_proportions
from .models import (
    EnvelopeRequest,
    EnvelopeResponse,
    FittedPeak,
    MixtureFitRequest,
    MixtureFitResponse,
    Peak,
    Proportion,
    ResidualPeak,
)

app = FastAPI(
    title="Isotope Envelope API",
    summary="同位素峰包络计算",
    description=(
        "以各元素最小质量数为基准，按独立抽样卷积生成分子的质量偏移分布；"
        "相同偏移合并概率，并按概率加权给出精确质量质心。"
    ),
    version="1.0.0",
)


@app.get("/")
def root():
    return {"service": "isotope-envelope", "endpoint": "POST /envelope", "docs": "/docs"}


@app.post("/envelope", response_model=EnvelopeResponse)
def envelope(request: EnvelopeRequest) -> EnvelopeResponse:
    elements = []
    for el in request.elements:
        base = min(iso.mass_number for iso in el.isotopes)
        isotopes = [
            (iso.mass_number - base, iso.exact_mass, iso.abundance)
            for iso in el.isotopes
        ]
        elements.append((el.count, isotopes))

    peaks, total = compute_envelope(elements)
    return EnvelopeResponse(
        peaks=[
            Peak(offset=offset, probability=str(prob), centroid=str(centroid))
            for offset, prob, centroid in peaks
        ],
        total_probability=str(total),
    )


@app.post("/mixture-fit", response_model=MixtureFitResponse)
def mixture_fit(request: MixtureFitRequest) -> MixtureFitResponse:
    # 各候选的绝对质量数归一化概率谱（复用包络计算；
    # 绝对质量数 = 各元素最小质量数 × 原子数之和 + 偏移）
    spectra = []
    for candidate in request.candidates:
        elements = []
        base = 0
        for el in candidate.elements:
            el_base = min(iso.mass_number for iso in el.isotopes)
            base += el_base * el.count
            isotopes = [
                (iso.mass_number - el_base, iso.exact_mass, iso.abundance)
                for iso in el.isotopes
            ]
            elements.append((el.count, isotopes))
        peaks, _ = compute_envelope(elements)
        spectrum: dict[int, float] = {}
        for offset, prob, _centroid in peaks:
            mass = base + offset
            spectrum[mass] = spectrum.get(mass, 0.0) + float(prob)
        total = sum(spectrum.values())
        spectra.append({mass: p / total for mass, p in spectrum.items()})

    # 实测归一化谱
    measured_total = sum(float(p.intensity) for p in request.measured)
    measured = {
        p.mass_number: float(p.intensity) / measured_total for p in request.measured
    }

    # 双方质量数并集（缺省补零）
    grid = sorted(set().union(*spectra, measured))
    index = {mass: i for i, mass in enumerate(grid)}
    columns = [[0.0] * len(grid) for _ in spectra]
    for j, spectrum in enumerate(spectra):
        for mass, prob in spectrum.items():
            columns[j][index[mass]] = prob
    b = [0.0] * len(grid)
    for mass, value in measured.items():
        b[index[mass]] = value

    proportions = fit_proportions(columns, b)

    fitted = [
        sum(columns[j][i] * proportions[j] for j in range(len(columns)))
        for i in range(len(grid))
    ]
    residuals = [b[i] - fitted[i] for i in range(len(grid))]
    rss = sum(r * r for r in residuals)

    return MixtureFitResponse(
        proportions=[
            Proportion(id=candidate.id, proportion=proportions[j])
            for j, candidate in enumerate(request.candidates)
        ],
        fitted=[
            FittedPeak(mass_number=mass, intensity=fitted[i])
            for i, mass in enumerate(grid)
        ],
        residuals=[
            ResidualPeak(mass_number=mass, residual=residuals[i])
            for i, mass in enumerate(grid)
        ],
        rss=rss,
    )
