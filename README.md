# Isotope Envelope API

同位素峰包络计算服务（Python + FastAPI）。以各元素最小质量数为基准，把每个原子看作
独立抽样，通过卷积生成分子的质量偏移分布；相同偏移合并概率，并按概率加权计算
精确质量质心。全程使用高精度 `Decimal`（80 位有效数字），不裁剪任何低概率中间项。

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

## 接口

- 计算接口：`POST http://localhost:8000/envelope`
- 交互式文档（Swagger UI）：`http://localhost:8000/docs`
- 服务信息：`GET http://localhost:8000/`

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
Dockerfile      镜像定义（启动时安装依赖）
docker-compose.yml
requirements.txt
```
