from opensuse/leap

WORKDIR /root
RUN zypper install -y git gcc autoconf automake libbz2-devel gsl-devel nasm yasm zlib-devel liblzma5 \
	xz-devel make gcc-c++

RUN zypper install -y libtool cmake which \
        java-17-openjdk \
        java-17-openjdk-devel \
	java-17-openjdk-headless \
	git-lfs \
	python311 \
	gnu_parallel

RUN cd /usr/bin/ && ln -s /usr/bin/python3.11 /usr/bin/python

#COPY build_deps.sh /root/build_deps.sh
#RUN ./build_deps.sh
