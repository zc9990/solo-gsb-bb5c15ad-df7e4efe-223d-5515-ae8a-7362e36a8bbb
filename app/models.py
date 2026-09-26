"""请求 / 响应模型与输入校验。

所有校验失败都会由 FastAPI 统一返回 422，并在响应体的
``detail[].loc`` 中给出出错的字段位置。
"""

from decimal import Decimal
from typing import Annotated, List, Optional

from pydantic import BaseModel, Field, field_validator, model_validator

MAX_MASS_NUMBER_SPAN = 4
MIN_TOTAL_ATOMS = 1
MAX_TOTAL_ATOMS = 60
MIN_CANDIDATES = 2
MAX_CANDIDATES = 8


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


class ElementsRequest(BaseModel):
    """含 ``elements`` 字段的请求公共部分（元素组成与总原子数校验）。"""

    elements: List[Element] = Field(..., min_length=1)

    @model_validator(mode="after")
    def check_total_atoms(self) -> "ElementsRequest":
        total = sum(el.count for el in self.elements)
        if not MIN_TOTAL_ATOMS <= total <= MAX_TOTAL_ATOMS:
            raise ValueError(
                f"总原子数必须在 {MIN_TOTAL_ATOMS}~{MAX_TOTAL_ATOMS} 之间，实际为 {total}"
            )
        return self


class EnvelopeRequest(ElementsRequest):
    """POST /envelope 请求体。"""


class MixtureCandidate(ElementsRequest):
    """一个候选分子：唯一 id + 元素组成（校验规则与 /envelope 相同）。"""

    id: str = Field(..., min_length=1, description="候选唯一标识")


class MixtureFitRequest(BaseModel):
    """POST /mixture-fit 请求体。"""

    candidates: List[MixtureCandidate] = Field(
        ...,
        min_length=MIN_CANDIDATES,
        max_length=MAX_CANDIDATES,
        description=f"{MIN_CANDIDATES}~{MAX_CANDIDATES} 个候选分子，id 唯一",
    )
    mass_number: List[Annotated[int, Field(ge=0)]] = Field(
        ..., min_length=1, description="实测谱整数质量数（非负、不重复）"
    )
    intensity: List[Annotated[Decimal, Field(ge=0)]] = Field(
        ..., min_length=1, description="实测谱强度（非负，总和大于 0）"
    )

    @model_validator(mode="after")
    def check_candidate_ids(self) -> "MixtureFitRequest":
        ids = [c.id for c in self.candidates]
        if len(set(ids)) != len(ids):
            raise ValueError("候选 id 存在重复值")
        return self

    @model_validator(mode="after")
    def check_measured(self) -> "MixtureFitRequest":
        if len(self.mass_number) != len(self.intensity):
            raise ValueError("mass_number 与 intensity 长度必须一致")
        if len(set(self.mass_number)) != len(self.mass_number):
            raise ValueError("实测 mass_number 存在重复值")
        total = sum(self.intensity, Decimal(0))
        if total <= 0:
            raise ValueError(f"实测强度总和必须大于 0，实际为 {total}")
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


class CandidateProportion(BaseModel):
    """一个候选的拟合比例。"""

    id: str
    proportion: float = Field(..., description="拟合比例（非负，全部比例之和为 1）")


class FittedIntensity(BaseModel):
    """拟合谱在一个质量数上的强度。"""

    mass_number: int
    intensity: float


class Residual(BaseModel):
    """一个质量数上「实测 − 拟合」的残差。"""

    mass_number: int
    residual: float


class MixtureFitResponse(BaseModel):
    """POST /mixture-fit 响应体。"""

    proportions: List[CandidateProportion] = Field(
        ..., description="按输入顺序的候选 id 与比例"
    )
    fitted: List[FittedIntensity] = Field(..., description="按质量数升序的拟合强度")
    residuals: List[Residual] = Field(
        ..., description="按质量数升序的残差（实测 − 拟合）"
    )
    sse: float = Field(..., description="残差平方和")
