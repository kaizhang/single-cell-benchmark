FROM python:3.7-slim
RUN apt-get update && apt-get install -y procps
RUN pip install --no-cache anndata scale==1.1.2 leidenalg