FROM rockylinux:9
RUN dnf -y install dnf-plugins-core epel-release
RUN dnf config-manager --set-enabled crb && dnf -y update
RUN dnf -y install python39 python39-wheel
RUN pip3 install --upgrade pip
RUN pip3 install scanpy
RUN pip3 install snapatac2