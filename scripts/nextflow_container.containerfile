from opensuse/leap


RUN zypper install -y \
	which \
	java-17-openjdk \
	java-17-openjdk-headless \
	podman

RUN curl -s https://get.nextflow.io | bash

