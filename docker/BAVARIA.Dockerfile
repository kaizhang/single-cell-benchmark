FROM rockylinux:8
RUN dnf -y install dnf-plugins-core epel-release
RUN dnf config-manager --set-enabled powertools && dnf -y update
RUN dnf -y install python38 python38-wheel
RUN pip3 install --upgrade pip setuptools
RUN pip3 install protobuf==3.20.*
RUN pip3 install https://github.com/BIMSBbioinfo/bavaria/archive/v0.1.0.zip