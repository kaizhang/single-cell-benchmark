FROM pytorch/pytorch:2.0.1-cuda11.7-cudnn8-runtime
RUN pip install --no-cache anndata numpy==1.23.* scale==1.1.2 leidenalg