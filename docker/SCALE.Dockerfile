FROM python:3.7-slim
RUN apt-get update && apt-get install -y --no-install-recommends gcc
RUN pip install torch torchvision torchaudio --extra-index-url https://download.pytorch.org/whl/cu113
RUN pip install anndata scale==1.1.2 leidenalg