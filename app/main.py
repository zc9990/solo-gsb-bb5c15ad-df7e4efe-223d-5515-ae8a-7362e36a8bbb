"""FastAPI 应用入口。"""

from fastapi import FastAPI

from .envelope import compute_envelope
from .models import EnvelopeRequest, EnvelopeResponse, Peak

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
