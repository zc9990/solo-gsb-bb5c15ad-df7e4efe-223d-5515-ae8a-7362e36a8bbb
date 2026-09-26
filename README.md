# Isotope Envelope API

同位素峰包络计算与混合谱拟合服务（Python + FastAPI）。

- `POST /envelope`：以各元素最小质量数为基准，把每个原子看作独立抽样，通过卷积生成
  分子的质量偏移分布；相同偏移合并概率，并按概率加权计算精确质量质心。全程使用高精度
  `Decimal`（80 位有效数字），不裁剪任何低概率中间项。
- `POST /mixture-fit`：给定 2~8 个候选分子与实测谱，求非负且和为 1 的候选比例，使
  拟合谱与实测谱的残差平方和全局最小；多解时取比例平方和最小者。拟合求解器为纯标准库
  实现（子集枚举 + KKT 最小范数解），**无新增第三方依赖**。

## 启动

单条命令（构建镜像并启动，容器启动时自动安装依赖）：

```bash
docker compose up --build
```

## 配置

| 配置项 | 默认值 | 说明 |
| --- | --- | --- |
| 端口映射 | `8000:8000` | 在 `docker-compose.yml` 的 `ports` 中修改，如 `"8080:8000"` 后访问 8080 端口 |

无其他必需配置（无需环境变量、数据库等）。

## 访问地址

服务启动后监听 `http://localhost:8000`（端口 8000）：

- 包络计算：`POST http://localhost:8000/envelope`
- 混合谱拟合：`POST http://localhost:8000/mixture-fit`
- 交互式文档（Swagger UI）：`http://localhost:8000/docs`
- 服务信息：`GET http://localhost:8000/`

## POST /envelope

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

## POST /mixture-fit

混合谱拟合：复用包络计算生成各候选的归一化概率谱，以**各元素最小质量数 × 原子数
之和 + 偏移**转换为绝对整数质量数；实测强度同样归一化，双方在质量数并集上补零。
然后求非负且和为 1 的候选比例，使拟合谱与实测谱的残差平方和（RSS）全局最小；
多解时取比例平方和最小者。数值误差不超过 1e-8。

### 请求体

```json
{
  "candidates": [
    {
      "id": "carbon",
      "elements": [
        {"name": "C", "count": 1, "isotopes": [
          {"mass_number": 12, "exact_mass": "12.0", "abundance": "0.9893"},
          {"mass_number": 13, "exact_mass": "13.003355", "abundance": "0.0107"}
        ]}
      ]
    },
    {
      "id": "oxygen",
      "elements": [
        {"name": "O", "count": 1, "isotopes": [
          {"mass_number": 16, "exact_mass": "15.994915", "abundance": "0.99757"},
          {"mass_number": 17, "exact_mass": "16.999131", "abundance": "0.00038"},
          {"mass_number": 18, "exact_mass": "17.999160", "abundance": "0.00205"}
        ]}
      ]
    }
  ],
  "measured": [
    {"mass_number": 12, "intensity": 2.9679},
    {"mass_number": 13, "intensity": 0.0321},
    {"mass_number": 16, "intensity": 6.98299},
    {"mass_number": 17, "intensity": 0.00266},
    {"mass_number": 18, "intensity": 0.01435}
  ]
}
```

字段说明：

- `candidates`：2～8 个候选分子。
  - `id`：候选唯一标识（字符串或整数），候选间不得重复；
  - `elements`：元素组成，模型与校验规则完全沿用 `/envelope`（含总原子数 1～60）。
- `measured`：实测谱，至少 1 个峰。
  - `mass_number`：整数质量数（≥ 0），不得重复；
  - `intensity`：强度（≥ 0），全部强度之和必须大于 0。

违反任一规则即返回 **422**（响应体 `detail[].loc` 标明字段位置）。

### 响应

```json
{
  "proportions": [
    {"id": "carbon", "proportion": 0.3},
    {"id": "oxygen", "proportion": 0.7}
  ],
  "fitted": [
    {"mass_number": 12, "intensity": 0.29679},
    {"mass_number": 13, "intensity": 0.00321},
    {"mass_number": 16, "intensity": 0.698299},
    {"mass_number": 17, "intensity": 0.000266},
    {"mass_number": 18, "intensity": 0.001435}
  ],
  "residuals": [
    {"mass_number": 12, "residual": 0.0},
    {"mass_number": 13, "residual": 0.0},
    {"mass_number": 16, "residual": 0.0},
    {"mass_number": 17, "residual": 0.0},
    {"mass_number": 18, "residual": 0.0}
  ],
  "rss": 1.23e-32
}
```

- `proportions`：按**输入顺序**的候选 `id` 与拟合比例（非负，全体和为 1）；
- `fitted`：按质量数升序的拟合强度（覆盖双方质量数并集，缺省为 0）；
- `residuals`：按质量数升序的残差（实测 − 拟合），与 `fitted` 同一质量数网格；
- `rss`：残差平方和（归一化谱尺度）。

### 示例调用

```bash
curl -X POST http://localhost:8000/mixture-fit \
  -H 'Content-Type: application/json' \
  -d '{"candidates":[{"id":"a","elements":[{"count":1,"isotopes":[{"mass_number":12,"exact_mass":"12.0","abundance":"0.9893"},{"mass_number":13,"exact_mass":"13.003355","abundance":"0.0107"}]}]},{"id":"b","elements":[{"count":1,"isotopes":[{"mass_number":16,"exact_mass":"15.994915","abundance":"1"}]}]}],"measured":[{"mass_number":12,"intensity":0.9893},{"mass_number":13,"intensity":0.0107},{"mass_number":16,"intensity":3.0}]}'
```

## 本地开发（不用 Docker）

```bash
pip install -r requirements.txt
uvicorn app.main:app --reload --port 8000
```

## 项目结构

```
app/
  main.py       FastAPI 入口与路由（/envelope、/mixture-fit）
  models.py     请求/响应模型与 422 校验
  envelope.py   卷积计算核心（高精度 Decimal）
  mixture.py    混合谱拟合求解器（纯标准库，无第三方依赖）
Dockerfile      镜像定义（启动时安装依赖）
docker-compose.yml
requirements.txt
```
