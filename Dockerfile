FROM python:3.12-slim

WORKDIR /app

# 仅拷贝代码与依赖清单；依赖在容器启动时安装（见 CMD）
COPY requirements.txt ./
COPY app ./app

EXPOSE 8000

# 启动时先安装依赖，再启动服务
CMD ["sh", "-c", "pip install --no-cache-dir -r requirements.txt && uvicorn app.main:app --host 0.0.0.0 --port 8000"]
