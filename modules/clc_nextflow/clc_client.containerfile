FROM opensuse/leap:latest

RUN zypper install -y gzip tar
ARG clc_servertools_install_sh=CLCServerCommandLineTools_23_0_4_64.sh

COPY $clc_servertools_install_sh /root/clcclient.sh
RUN chmod u+x /root/clcclient.sh
# install clcserver without confirmation prompts
RUN /root/clcclient.sh -q -dir /opt/clcclient

CMD ["/opt/clcclient/clcserver"]
