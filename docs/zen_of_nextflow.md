# Conventions for good Nextflow code

these are my conventions for writing good nextflow code.

they are loosely inspired by the zen of python and build on my practical experience from writing nextflow pipelines.

those apply to code, documentation and data.

## In order of priority, with highest priority first:

* simplicity > performance
* concise > verbosity
* dont use code styling (ascii-art)
* good style comes through form follows function
* keep high coupled interactions local (code and data-wise)
* dont repeat yourself (DRY)
* but repeat yourself if it avoids alot of complexity
* optimize late
* when optimizing, optimize separately and at a lower level
* explicit > implicit
* flat > nested
* fail fast and hard
* check inputs where the inputs are used
* if you need to check inputs earlier, scan the input checks and generate a checker (in other words, dont use hardcoded checks as late as necessary)
* therefore fail where the error occurs, as this improves debugging alot
* functional, immutable, stateless > inplace stateful
* use types for correctness > unittesting first > runtime validation
* hackable code > extendable code > plugins > configuration
