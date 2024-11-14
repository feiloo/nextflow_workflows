from pipeline_task

RUN zypper install -y \
	which \
	java-17-openjdk \
	java-17-openjdk-headless


RUN curl -s https://get.nextflow.io | bash
RUN cp nextflow /usr/local/bin/

RUN zypper install -y meson

#RUN git clone --recursive https://github.com/feiloo/nextflow_workflows.git
RUN mkdir -p /usr/local/lib/nextflow
COPY . /root/nextflow_workflows
RUN cd nextflow_workflows && meson setup ../nextflow_workflows_build && cd ../nextflow_workflows_build && meson install

ENV NEXTFLOW_MODULES /usr/local/lib/nextflow/nextflow_workflows/modules
ENV NGS_REFERENCE_DIR /usr/share/nextflow_pipeline_reference
ENV NEXTFLOW_STOREDIR /var/cache/nextflow_pipeline_storedir
ENV NEXTFLOW_OUTPUTDIR_CUSTOM /root/nextflow_outputdir
ENV NEXTFLOW_WORKDIR_CUSTOM /root/nextflow_workdir
ENV NAS_IMPORT_DIR ''
ENV NAS_EXPORT_DIR ''
