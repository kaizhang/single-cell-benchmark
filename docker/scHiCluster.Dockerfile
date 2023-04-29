FROM python:3.10-slim
RUN apt-get update && apt-get install -y procps r-base
RUN pip install --no-cache anndata opencv-python schicluster==1.3.2