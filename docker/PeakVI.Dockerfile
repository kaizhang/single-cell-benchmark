FROM pytorch/pytorch:1.13.1-cuda11.6-cudnn8-runtime
RUN pip install --no-cache-dir anndata scvi-tools==0.19.0