# Isotope Envelope API

同位素峰包络计算与混合谱拟合服务（Python + FastAPI）。以各元素最小质量数为基准，把每个原子看作
独立抽样，通过卷积生成分子的质量偏移分布；相同偏移合并概率，并按概率加权计算
精确质量质心。全程使用高精度 `Decimal`（80 位有效数字），不裁剪任何低概率中间项。

- `POST /envelope`：计算单个分子组成的同位素峰包络。
- `POST /mixture-fit`：给定 2~8 个候选分子与一条实测谱，求非负且和为 1 的候选比例，
  使拟合谱与实测谱的残差平方和全局最小；多解时取比例平方和最小者。

## 启动

单条命令（构建镜像并启动，依赖在镜像构建阶段安装完成）：

```bash
docker compose up --build
```

启动后访问地址：http://localhost:8000

## 配置

| 配置项 | 默认值 | 说明 |
| --- | --- | --- |
| 端口映射 | `8000:8000` | 在 `docker-compose.yml` 的 `ports` 中修改，如 `"8080:8000"` 后访问 http://localhost:8080 |

无其他必需配置（无需环境变量、数据库等）。

## 接口

- 包络计算：`POST http://localhost:8000/envelope`
- 混合谱拟合：`POST http://localhost:8000/mixture-fit`
- 交互式文档（Swagger UI）：`http://localhost:8000/docs`
- 服务信息：`GET http://localhost:8000/`

## 包络计算 `POST /envelope`

### 请求体

```json
{
  "elements": [
    {
      "name": "H",
      "count": 2,
      "isotopes": [
        {"mass_number": 1, "exact_mass": "1.007825", "abundance": "0.999885"},
        {"mass_number": 2, "exact_mass": "2.014102", "abundance": "0.000115"}
      ]
    },
    {
      "name": "O",
      "count": 1,
      "isotopes": [
        {"mass_number": 16, "exact_mass": "15.994915", "abundance": "0.99757"},
        {"mass_number": 17, "exact_mass": "16.999131", "abundance": "0.00038"},
        {"mass_number": 18, "exact_mass": "17.999160", "abundance": "0.00205"}
      ]
    }
  ]
}
```

字段说明：

- `elements[].count`：该元素的原子数，非负整数；`name` 为可选标识。
- `elements[].isotopes[].mass_number`：整数质量数（≥ 0）。
- `elements[].isotopes[].exact_mass`：十进制精确质量（> 0），建议用字符串以保证精度，JSON 数字也可。
- `elements[].isotopes[].abundance`：丰度（≥ 0），同一元素的丰度合计必须恰好为 1。

校验规则（违反即返回 **422**，响应体 `detail[].loc` 标明字段位置）：

- 总原子数（所有元素 `count` 之和）必须在 **1～60** 之间；
- 同一元素的质量数跨度（最大 − 最小）不得超过 **4**，且质量数不得重复；
- 丰度非负且逐元素合计为 1；原子数为非负整数；同位素列表至少 1 项。

### 响应

```json
{
  "peaks": [
    {"offset": 0, "probability": "0.99734057209286325", "centroid": "18.010565"},
    {"offset": 1, "probability": "0.00060932731929900", "centroid": "19.015556976574727755648448779072572347699547126390759133530927437498092377814543"}
  ],
  "total_probability": "1.00000000000000000"
}
```

- `peaks`：按质量偏移 `offset` 升序排列，仅包含概率非零的峰；
  - `probability`：该偏移的合并概率（十进制字符串，完整保留精度）；
  - `centroid`：该偏移下所有同位素组合按概率加权的精确质量质心。
- `total_probability`：全部峰的概率总和（理论值为 1，返回实际计算值）。

### 示例调用

```bash
curl -X POST http://localhost:8000/envelope \
  -H 'Content-Type: application/json' \
  -d '{"elements":[{"count":1,"isotopes":[{"mass_number":12,"exact_mass":"12.0","abundance":"0.9893"},{"mass_number":13,"exact_mass":"13.003355","abundance":"0.0107"}]}]}'
```

## 混合谱拟合 `POST /mixture-fit`

### 请求体

```json
{
  "candidates": [
    {
      "id": "C",
      "elements": [
        {
          "count": 1,
          "isotopes": [
            {"mass_number": 12, "exact_mass": "12.0", "abundance": "0.99"},
            {"mass_number": 13, "exact_mass": "13.003355", "abundance": "0.01"}
          ]
        }
      ]
    },
    {
      "id": "N",
      "elements": [
        {
          "count": 1,
          "isotopes": [
            {"mass_number": 14, "exact_mass": "14.003074", "abundance": "0.996"},
            {"mass_number": 15, "exact_mass": "15.000109", "abundance": "0.004"}
          ]
        }
      ]
    }
  ],
  "mass_number": [12, 13, 14, 15],
  "intensity": ["495", "5", "498", "2"]
}
```

字段说明：

- `candidates`：**2~8** 个候选分子；每个候选含：
  - `id`：非空字符串，候选间必须唯一；
  - `elements`：元素组成，结构与校验规则与 `/envelope` 完全一致。
- `mass_number`：实测谱整数质量数列表，逐项对应 `intensity`，每个值必须 ≥ 0 且不重复。
- `intensity`：实测强度列表，每个值必须 ≥ 0（建议用字符串保证十进制精度），且总和必须 > 0。

以上任何校验失败均返回 **422**，响应体 `detail[].loc` 标明字段位置。

### 计算口径

1. 复用 `/envelope` 的包络计算得到每个候选的质量偏移分布，以
   **各元素最小质量数 × 原子数之和 + 偏移** 转换为绝对整数质量数谱；
2. 各候选谱按概率总和归一化，实测谱按强度总和归一化；
3. 在候选与实测所有质量数的**并集**网格上，缺失的一侧补零对齐；
4. 求比例向量 `x`（`xᵢ ≥ 0`、`Σxᵢ = 1`），使 `Σ(拟合强度 − 实测强度)²` 全局最小
   （凸二次规划，枚举支撑集解 KKT 条件得到全局最优）；最优解不唯一时，取 `Σxᵢ²` 最小者。

### 响应

```json
{
  "proportions": [
    {"id": "C", "proportion": 0.5},
    {"id": "N", "proportion": 0.5}
  ],
  "fitted": [
    {"mass_number": 12, "intensity": 0.495},
    {"mass_number": 13, "intensity": 0.005},
    {"mass_number": 14, "intensity": 0.498},
    {"mass_number": 15, "intensity": 0.002}
  ],
  "residuals": [
    {"mass_number": 12, "residual": 0.0},
    {"mass_number": 13, "residual": 0.0},
    {"mass_number": 14, "residual": 0.0},
    {"mass_number": 15, "residual": 0.0}
  ],
  "sse": 0.0
}
```

- `proportions`：按**输入顺序**给出的候选 id 与比例；非负、总和为 1。
- `fitted`：按质量数**升序**给出的拟合强度（与实测同口径，均为归一化后的值）。
- `residuals`：同序的残差（**实测 − 拟合**）。
- `sse`：残差平方和。数值结果的误差不超过 1e-8。

### 示例调用

```bash
curl -X POST http://localhost:8000/mixture-fit \
  -H 'Content-Type: application/json' \
  -d '{"candidates":[{"id":"C","elements":[{"count":1,"isotopes":[{"mass_number":12,"exact_mass":"12.0","abundance":"0.99"},{"mass_number":13,"exact_mass":"13.003355","abundance":"0.01"}]}]},{"id":"N","elements":[{"count":1,"isotopes":[{"mass_number":14,"exact_mass":"14.003074","abundance":"0.996"},{"mass_number":15,"exact_mass":"15.000109","abundance":"0.004"}]}]}],"mass_number":[12,13,14,15],"intensity":["495","5","498","2"]}'
```

## 本地开发（不用 Docker）

```bash
pip install -r requirements.txt
uvicorn app.main:app --reload --port 8000
```

## 项目结构

```
app/
  main.py       FastAPI 入口与路由
  models.py     请求/响应模型与 422 校验
  envelope.py   卷积计算核心（高精度 Decimal）
  mixture.py    混合谱拟合（归一化对齐 + 单纯形约束凸二次规划）
Dockerfile      镜像定义（构建阶段安装依赖）
docker-compose.yml
requirements.txt
```
