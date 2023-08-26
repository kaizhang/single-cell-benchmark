FROM pytorch/pytorch:2.0.1-cuda11.7-cudnn8-runtime
RUN pip install --no-cache-dir anndata scvi-tools==1.0.3