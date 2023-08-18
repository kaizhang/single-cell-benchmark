FROM python:3.10-slim
RUN apt-get update && apt-get install -y procps
RUN pip install --no-cache scanpy harmonypy==0.0.9 snapatac2==2.3.1