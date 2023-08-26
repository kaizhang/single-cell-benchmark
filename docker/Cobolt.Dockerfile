FROM pytorch/pytorch:2.0.1-cuda11.7-cudnn8-runtime
RUN apt-get update && apt-get install -y git
RUN pip install --no-cache anndata
RUN pip install --no-cache git+https://github.com/epurdom/cobolt.git@70b6ff1365c4fbd2161e4297f8455455e9e87354