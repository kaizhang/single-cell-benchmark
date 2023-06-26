FROM python:3.10-slim
RUN apt-get update && apt-get install -y procps
RUN pip install --no-cache scanpy snapatac2==2.3.1