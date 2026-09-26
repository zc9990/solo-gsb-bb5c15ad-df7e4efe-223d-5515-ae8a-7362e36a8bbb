"""请求 / 响应模型与输入校验。

所有校验失败都会由 FastAPI 统一返回 422，并在响应体的
``detail[].loc`` 中给出出错的字段位置。
"""

from decimal import Decimal
from typing import List, Optional, Union

from pydantic import BaseModel, Field, field_validator, model_validator

MAX_MASS_NUMBER_SPAN = 4
MIN_TOTAL_ATOMS = 1
MAX_TOTAL_ATOMS = 60


class Isotope(BaseModel):
    """一种同位素：整数质量数、十进制精确质量、丰度。"""

    mass_number: int = Field(..., ge=0, description="整数质量数")
    exact_mass: Decimal = Field(..., gt=0, description="十进制精确质量")
    abundance: Decimal = Field(..., ge=0, description="丰度（非负）")


class Element(BaseModel):
    """一种元素：非负整数原子数 + 同位素列表。"""

    name: Optional[str] = Field(None, description="元素符号（可选，仅作标识）")
    count: int = Field(..., ge=0, description="原子数（非负整数）")
    isotopes: List[Isotope] = Field(..., min_length=1)

    @field_validator("isotopes")
    @classmethod
    def check_mass_numbers(cls, isotopes: List[Isotope]) -> List[Isotope]:
        numbers = [iso.mass_number for iso in isotopes]
        if len(set(numbers)) != len(numbers):
            raise ValueError("mass_number 存在重复值")
        if max(numbers) - min(numbers) > MAX_MASS_NUMBER_SPAN:
            raise ValueError(
                f"同一元素质量数跨度不能超过 {MAX_MASS_NUMBER_SPAN}"
            )
        return isotopes

    @model_validator(mode="after")
    def check_abundance_sum(self) -> "Element":
        total = sum((iso.abundance for iso in self.isotopes), Decimal(0))
        if total != 1:
            raise ValueError(f"同一元素的丰度合计必须为 1，实际为 {total}")
        return self


class EnvelopeRequest(BaseModel):
    """POST /envelope 请求体。"""

    elements: List[Element] = Field(..., min_length=1)

    @model_validator(mode="after")
    def check_total_atoms(self) -> "EnvelopeRequest":
        total = sum(el.count for el in self.elements)
        if not MIN_TOTAL_ATOMS <= total <= MAX_TOTAL_ATOMS:
            raise ValueError(
                f"总原子数必须在 {MIN_TOTAL_ATOMS}~{MAX_TOTAL_ATOMS} 之间，实际为 {total}"
            )
        return self


class Peak(BaseModel):
    """一个峰：质量偏移、非零概率、概率加权精确质量质心。

    概率与质心以十进制字符串输出，完整保留精度。
    """

    offset: int
    probability: str
    centroid: str


class EnvelopeResponse(BaseModel):
    """POST /envelope 响应体。"""

    peaks: List[Peak] = Field(..., description="按偏移升序的非零概率峰")
    total_probability: str = Field(..., description="全部峰的概率总和")


class Candidate(EnvelopeRequest):
    """一个候选分子：唯一 id + 元素组成。

    元素模型与校验（原子数、质量数跨度、丰度合计等）完全沿用 /envelope。
    """

    id: Union[str, int] = Field(..., description="候选唯一标识（字符串或整数）")


class MeasuredPeak(BaseModel):
    """实测谱的一个峰：非负整数质量数 + 非负强度。"""

    mass_number: int = Field(..., ge=0, description="整数质量数（非负）")
    intensity: Decimal = Field(..., ge=0, description="强度（非负）")


class MixtureFitRequest(BaseModel):
    """POST /mixture-fit 请求体：2~8 个候选 + 实测谱。"""

    candidates: List[Candidate] = Field(..., min_length=2, max_length=8)
    measured: List[MeasuredPeak] = Field(..., min_length=1)

    @model_validator(mode="after")
    def check_unique_ids(self) -> "MixtureFitRequest":
        ids = [c.id for c in self.candidates]
        if len(set(ids)) != len(ids):
            raise ValueError("候选 id 存在重复")
        return self

    @model_validator(mode="after")
    def check_measured(self) -> "MixtureFitRequest":
        numbers = [p.mass_number for p in self.measured]
        if len(set(numbers)) != len(numbers):
            raise ValueError("实测 mass_number 存在重复值")
        total = sum((p.intensity for p in self.measured), Decimal(0))
        if total <= 0:
            raise ValueError(f"实测强度总和必须大于 0，实际为 {total}")
        return self


class Proportion(BaseModel):
    """一个候选的拟合比例。"""

    id: Union[str, int] = Field(..., description="候选 id（与输入一致）")
    proportion: float = Field(..., description="拟合比例（非负，全体和为 1）")


class FittedPeak(BaseModel):
    """拟合谱的一个点：质量数 + 拟合强度。"""

    mass_number: int
    intensity: float


class ResidualPeak(BaseModel):
    """一个质量数上的残差（实测 − 拟合）。"""

    mass_number: int
    residual: float


class MixtureFitResponse(BaseModel):
    """POST /mixture-fit 响应体。"""

    proportions: List[Proportion] = Field(..., description="按输入顺序的候选 id 与比例")
    fitted: List[FittedPeak] = Field(..., description="按质量数升序的拟合强度")
    residuals: List[ResidualPeak] = Field(
        ..., description="按质量数升序的残差（实测 − 拟合）"
    )
    rss: float = Field(..., description="残差平方和")
