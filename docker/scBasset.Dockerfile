FROM tensorflow/tensorflow:2.11.0-gpu
RUN apt-get update && apt-get install -y procps git libbz2-dev liblzma-dev
RUN git clone https://github.com/calico/scBasset.git
RUN cd scBasset && pip install --no-cache .
RUN cd .. && rm -rf scBasset