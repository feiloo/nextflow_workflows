FROM opensuse/leap

RUN zypper install -y gzip tar

COPY CLCServerCommandLineTools_23_0_4_64.sh /root/clcclient.sh

# install clcserver iwthout confirmation prompts
RUN /root/clcclient.sh -q -dir /opt/clcclient

CMD ["/opt/clcclient/clcserver"]
