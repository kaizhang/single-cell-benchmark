FROM python:3.10-slim
RUN apt-get update && apt-get install -y procps
RUN pip install --no-cache scanpy muon mofapy2==0.7.0