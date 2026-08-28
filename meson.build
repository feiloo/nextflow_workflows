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

dwgsim_src = meson.project_source_root() / 'third_party/DWGSIM/.'
dwgsim_build_dir = meson.current_build_dir() / 'third_party/DWGSIM'

build_dwgsim = run_target('build_dwgsim',
    command: ['bash', '-c', 
      '''set -euo pipefail 
      mkdir -p third_party
      cp -r @0@ @1@ 
      cd third_party/DWGSIM 
      sed -i "s/-lcurses/-lncurses/g" samtools/Makefile 
      make -j $(nproc)
      '''.format(dwgsim_src, dwgsim_build_dir)
      ],
)

dwgsim_bin = custom_target('dwgsim',
    output: 'dwgsim',   # will be placed in dwgsim_build_dir
    build_by_default: false,
    command: [
        'bash', '-c',
        '''
        set -euo pipefail
        mkdir -p third_party
        # Copy entire submodule directory into build dir
        cp -r "@0@" "@1@"   # @0 = source, @1 = destination (dwgsim_build_dir)
        cd "@1@"
        # Apply patch
        sed -i 's/-lcurses/-lncurses/' samtools/Makefile
        # Build
        make -j $(nproc)
        '''
        .format(dwgsim_src, dwgsim_build_dir)
    ],
    install: false
)

# Run target – calls the bash script with test_work as both the working dir and output dir
# this generates the expected testdata (stub) in meson_builddir/test_work
run_target('generate-wgs',
  command: [
    'mkdir','-p', meson.current_build_dir() / 'scripts','&&',
    'cp', '-r', meson.project_source_root() / 'scripts/generate_stub_data.sh', meson.current_build_dir() / 'scripts', '&&',
    './scripts/generate_stub_data.sh',
    ref_genome,
    test_work,       # output directory
    output_prefix,
    test_work        # persistent temporary directory (same as output dir)
  ],
  depends: [dwgsim_bin]
)
