this lists some issues that became apparent when processing wgs samples.

this is not an exhaustive list

one fastq sample was 800GB gzip compressed per read, 3.8TB uncompressed

`seqtk sample` (https://github.com/lh3/seqtk) oomed on with above 700GB ram
`seqtk sample -2` with its flag for two passes for slower but lower memory execution, didnt finish within 9 days

fastp as in this pipeline, oomed with above 800GB ram

the sample didnt go further than ooming on fastp in this pipeline


gzip decoding can half-silently fail by `ignoring trailing garbage` on truncated or otherwise corrupted files.

this can lead to loosing data if not noticed.


gzip -l, showing uncompressed size and compressionratio, seems to output nonsense numbers in some cases, like the above


compressed filesize is a bad/worse heuristic for required compute time and memory, than expected.
compression ratio can vary more than expected.

bgzip is significantly faster, even with 1 thread than gzip.
