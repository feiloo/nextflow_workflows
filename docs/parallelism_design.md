### our paralellism strategy:

*subject to change*

use nextflow processes (dataflow) to use multinode+multiprocessing parallelise concurrently across samples/files

we work to enable multinode, but for now we just use 1 node

for now, we avoid splitting files for more parallelism

instead we use tool-level threading

so we run multiple nf-processes

each nf-process runs a single tool process 

each tool process runs the tool with the tools parallelism option


with the current filesize to count ratio: 1GB to 5 GB per file and about 5-50 files per workflow on ~100 cores this works pretty well

and most importantly it is quite simple to reason about
