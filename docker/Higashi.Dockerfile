FROM pytorch/pytorch:1.13.1-cuda11.6-cudnn8-runtime
RUN apt-get update && apt-get install -y procps git libbz2-dev liblzma-dev
RUN git clone https://github.com/ma-compbio/Higashi/ && cd Higashi && python setup.py install
RUN rm -r Higashi