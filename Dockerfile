FROM python:3.12-slim

WORKDIR /app

# 先拷贝依赖清单并安装（构建阶段完成，docker compose up --build 即可就绪）
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# 再拷贝应用代码
COPY app ./app

EXPOSE 8000

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]
