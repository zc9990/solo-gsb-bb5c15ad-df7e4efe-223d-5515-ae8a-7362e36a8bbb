"""请求 / 响应模型与输入校验。

所有校验失败都会由 FastAPI 统一返回 422，并在响应体的
``detail[].loc`` 中给出出错的字段位置。
"""

from decimal import Decimal
from typing import List, Optional

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
