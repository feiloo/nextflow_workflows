# clc nextflow wrapper module

at this time, this wrapper only calls a peculiar one of our internal clc workflows.

this module calls a clc workflow from within a nextflow workflow.
therefore it allows combining clc workflows with nextflow ones.

to use this, you need a clcserver setup and clc workflows


to add more clc workflows, one can simply add nextflow process definitions for each one, like a library.


## building the clc container
you need to download the clc server commandline tools to [the clc website](https://digitalinsights.qiagen.com/products-overview/discovery-insights-portfolio/enterprise-ngs-solutions/clc-server-command-line-tools/)


then you can build the container, see `build_container.sh`

now you can test (-stub) or run the workflow like this:

```
nextflow run $NEXTFLOW_MODULES/clc_nextflow -c $NEXFLOW_CONFIGS/your_config.config --samplesheet path_to_your_samplesheet.csv -resume -stub
```
