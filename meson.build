project('nextflow_workflows')

subdir('modules')

# generate some stub wgs data:

ref_dir = get_option('reference_dir')
if ref_dir == ''
  error('ref_genome option must be set to a valid FASTA file')
endif

ref_genome = ref_dir / 'Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta'


# Working directory inside the build tree (both temporary and output go here)
test_work = meson.current_build_dir() / 'test_work'

# Output prefix
output_prefix = 'testdata'

# Run target – calls the bash script with test_work as both the working dir and output dir
run_target('generate-wgs',
  command: [
    files('scripts/generate_stub_data.sh'),
    ref_genome,
    test_work,       # output directory
    output_prefix,
    test_work        # persistent temporary directory (same as output dir)
  ],
  depends: []
)
