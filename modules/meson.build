# nextflow_bin = find_program('nextflow')
nextflow_bin = find_program('nextflow-25.04.7')
nextflow_libdir = '/usr/local/lib/nextflow'
test_datadir = get_option('test_datadir')
enable_proprietary = get_option('enable_proprietary')



subdir('fusions')

nextpipe_sources = [
'nextpipe/resources/usr/bin/transcript_list_08_08_2024.py',
'nextpipe/resources/usr/bin/organize_variant_table.py',
'nextpipe/resources/usr/bin/variantlist_utils.py',
]

if enable_proprietary
  clc_workflow_name = get_option('clc_workflow_name')
  clc_import_dir = get_option('clc_import_dir')
  clc_export_dir = get_option('clc_export_dir')
  clc_destdir = get_option('clc_destdir')
  subdir('clc_nextflow')
endif

subdir('sequence_alignment')
subdir('biomarker')

sources = ['oncoscanner'/'user.config']
testdata = []

create_nextpipe_dir = custom_target('create_nextpipe_dir',
	command : ['mkdir', '-p', meson.current_build_dir() / 'nextpipe'/'bin', 
		'&&', 
		'touch', 'nextpipe_testfile'
		],
	output : ['nextpipe_testfile']
	)

copy_nextpipe_bin = custom_target('copy_nextpipe_bin',
	input : nextpipe_sources,
	output : ['organize_variant_table.py','variantlist_utils.py', 'transcript_list_08_08_2024'],
	command : ['mkdir', '-p', '@OUTDIR@/nextpipe/resources/usr/bin', '&&', 'cp', '-r', '@INPUT@', '@OUTDIR@/nextpipe/resources/usr/bin/'],
	install : false,
	depends : [create_nextpipe_dir]
	)

copy_nextpipe_sources = custom_target('copy_nextpipe_main',
	input : 'nextpipe/main.nf',
	output : 'main.nf',
	command : ['cp', '-r', '@INPUT@', '@OUTDIR@/nextpipe/'],
	install : false,
	depends : [copy_nextpipe_bin]
	)


copy_main_sources = custom_target('copy_main_sources',
	input : sources + testdata,
	output : 'user.config',
	command : ['cp', '-r', '@INPUT@', '@OUTDIR@/oncoscanner'],
	install : false,
        depends : [biomarker]
	)


# collected by running "find" for .nf environment.yml meta.yml .sh .py and .config files
# TODO: add all used nf source files for proper buildsystem detection
variantinterpretation_files = [
  'variantinterpretation_vendor/variantinterpretation.nf',
  ]

# just builds/copies the variantinterpretation, since it doesnt have meson
variantinterpretation = custom_target('variantinterpretation',
	input : variantinterpretation_files,
	output : 'variantinterpretation_vendor',
	command : ['cp', '-rT', '@CURRENT_SOURCE_DIR@/variantinterpretation_vendor', '@OUTPUT@/'],
	)

# TODO: think about where assets would be better and store and copy from there
assets_files = [
  '../assets/QIAseqExomehg380720.bed',
]

# just builds/copies the variantinterpretation, since it doesnt have meson
assets = custom_target('assets',
	input : assets_files,
	output : 'QIAseqExomehg380720.bed',
	command : ['cp', '@INPUT@', '@OUTDIR@/variantinterpretation_vendor/@PLAINNAME@'],
	)

# need to include here, so that variantinterpretation is defined
subdir('oncoscanner')

env_vars = ['NEXTFLOW_MODULES='+meson.current_build_dir(), 
              'NAS_IMPORT_DIR='+meson.current_build_dir()/'nas_import', 
              'NAS_EXPORT_DIR='+meson.current_build_dir()/'nas_export']

test('nextflow_workflow_variantinterpretation',
	nextflow_bin,
	is_parallel : false,
	suite : ['nextflow', 'workflow'],
	env: env_vars,
	args : [
	'run',
	meson.current_build_dir()+'/oncoscanner/main.nf',
	'-c', meson.current_build_dir()/ 'oncoscanner' / 'user.config',
	'--workflow_variation', 'variantinterpretation',
	'--samplesheet', test_datadir / 'test_variantinterpretation_samplesheet.csv',
	],
	depends: [copy_main_sources, variantinterpretation, oncoscanner],
	timeout: 6000
	)

test('nextflow_workflow_fusions',
	nextflow_bin,
	is_parallel : false,
	suite : ['nextflow', 'workflow'],
	env: env_vars,
	args : [
	'run',
	meson.current_build_dir()+'/oncoscanner/main.nf',
	'-c', meson.current_build_dir()/ 'oncoscanner' / 'user.config',
	'--workflow_variation', 'fusions',
	'--samplesheet', test_datadir / 'test_fusions_samplesheet.csv',
	],
	depends: [copy_main_sources, fusions, oncoscanner],
	timeout: 600
	)


test('nextflow_workflow_sequence_alignment',
	nextflow_bin,
	is_parallel : false,
	suite : ['nextflow', 'workflow'],
	env: env_vars,
	args : [
	'run',
	meson.current_build_dir()+'/oncoscanner/main.nf',
	'-c', meson.current_build_dir()/'oncoscanner' / 'user.config',
	'-profile', 'standard',
	'--workflow_variation', 'sequence_alignment',
	'--samplesheet', test_datadir/'test_sequence_alignment_samplesheet.csv',
        '--samplesheet_flavor','modal',
        '--hash_db', test_datadir/'md5sum.txt',
	],
	depends: [copy_main_sources, variantinterpretation, sequence_alignment_nextflow, oncoscanner],
	timeout: 360000
	)

test('nextflow_workflow_sequence_alignment_smoke',
	nextflow_bin,
	is_parallel : false,
	suite : ['nextflow', 'workflow', 'smoke'],
	env: env_vars,
	args : [
	'run',
	meson.current_build_dir()+'/oncoscanner/main.nf',
	'-c', meson.current_build_dir()/'oncoscanner' / 'user.config',
	'-profile', 'standard',
	'--workflow_variation', 'sequence_alignment',
        # '--samplesheet', meson.current_build_dir() / '..' / 'test_work' / 'samplesheet.csv',
        '--samplesheet', 'test_work' / 'samplesheet.csv',
        '--samplesheet_flavor','modal',
        '--input_dir', meson.current_build_dir() / '..' / 'test_work', 
        '--hash_db', meson.current_build_dir() / '..' / 'test_work' / 'md5sum.txt',
	],
	depends: [copy_main_sources, variantinterpretation, sequence_alignment_nextflow, oncoscanner, assets],
	timeout: 360000
	)

test('nextflow_workflow_sequence_alignment_local',
	nextflow_bin,
	is_parallel : false,
	suite : ['nextflow', 'workflow'],
	env: env_vars,
	args : [
	'run',
	meson.current_build_dir()+'/oncoscanner/main.nf',
	'-c', meson.current_build_dir()/'oncoscanner' / 'user.config',
	'-profile', 'local_nas_workdir',
	'--workflow_variation', 'sequence_alignment',
	'--samplesheet', test_datadir/'test_sequence_alignment_samplesheet.csv',
        '--samplesheet_flavor','modal',
        '--hash_db', test_datadir/'md5sum.txt',
	],
	depends: [copy_main_sources, variantinterpretation, sequence_alignment_nextflow, oncoscanner],
	timeout: 360000
	)

test('nextflow_workflow_sequence_alignment_hpc',
	nextflow_bin,
	is_parallel : false,
	suite : ['nextflow', 'workflow'],
	env: env_vars,
	args : [
	'run',
	meson.current_build_dir()+'/oncoscanner/main.nf',
	'-c', meson.current_build_dir()/'oncoscanner' / 'user.config',
	'-profile', 'hpc',
	'--workflow_variation', 'sequence_alignment',
	'--samplesheet', test_datadir/'test_sequence_alignment_samplesheet.csv',
        '--samplesheet_flavor','modal',
        '--hash_db', test_datadir/'md5sum.txt',
	],
	depends: [copy_main_sources, variantinterpretation, sequence_alignment_nextflow, oncoscanner],
	timeout: 360000
	)

test('nextflow_workflow_sequence_align_interpret',
	nextflow_bin,
	is_parallel : false,
	suite : ['nextflow', 'workflow'],
	env: env_vars,
	args : [
	'run',
	meson.current_build_dir()+'/oncoscanner/main.nf',
	'-c', meson.current_build_dir()/'oncoscanner' / 'user.config',
	'-profile', 'standard',
	'--workflow_variation', 'align_interpret',
        '--library_type', 'wgs',
	'--samplesheet', test_datadir/'test_sequence_alignment_samplesheet.csv',
        '--samplesheet_flavor','modal',
        '--input_dir', test_datadir,
        '--hash_db', test_datadir/'md5sum.txt',
	],
	depends: [copy_main_sources, variantinterpretation, sequence_alignment_nextflow, biomarker, oncoscanner],
	timeout: 360000
	)

test('nextflow_workflow_sequence_align_interpret_smoke',
	nextflow_bin,
	is_parallel : false,
	suite : ['nextflow', 'workflow', 'smoke'],
	env: env_vars,
	args : [
	'run',
	meson.current_build_dir()+'/oncoscanner/main.nf',
	'-c', meson.current_build_dir()/'oncoscanner' / 'user.config',
	'-profile', 'standard',
	'--workflow_variation', 'align_interpret',
        '--library_type','wgs',
        # '--samplesheet', meson.current_build_dir() / '..' / 'test_work' / 'samplesheet.csv',
        '--samplesheet', 'test_work' / 'samplesheet.csv',
        '--samplesheet_flavor','modal',
        '--input_dir', meson.current_build_dir() / '..' / 'test_work', 
        '--hash_db', meson.current_build_dir() / '..' / 'test_work' / 'md5sum.txt',
        '--bedfile', meson.current_build_dir() / 'variantinterpretation_vendor/QIAseqExomehg380720.bed',
	],
	depends: [copy_main_sources, variantinterpretation, sequence_alignment_nextflow, biomarker, oncoscanner, assets],
	timeout: 360000
	)



if enable_proprietary
  create_testdirs = custom_target('create_testdirs',
          command : ['mkdir', '-p', meson.current_build_dir() / 'nas_import', 
                  meson.current_build_dir() / 'nas_export', 
                  '&&', 'touch', 'nas_testfile'
                  ],
          output : ['nas_testfile']
          )

  test('nextflow_workflow_clc',
          nextflow_bin,
          is_parallel : false,
          suite : ['nextflow', 'workflow'],
          env : env_vars,
          args : [
          'run',
          meson.current_build_dir()+'/oncoscanner/main.nf',
          '-c', meson.current_build_dir() / 'oncoscanner' / 'user.config',
          '--workflow_variation', 'clc',
          '--samplesheet', test_datadir/'test_clc_samplesheet.csv',
          '-with-conda',
          '-stub'
          ],
          depends: [clc_nextflow, oncoscanner, create_testdirs, copy_sources, copy_nextpipe_sources, sequence_alignment_nextflow, variantinterpretation],
          timeout: 6000
          )

  test('nextflow_workflow_clc_dna_only_full',
          nextflow_bin,
          is_parallel : false,
          suite : ['nextflow', 'workflow', 'full'],
          env : env_vars,
          args : [
          'run',
          meson.current_build_dir()+'/oncoscanner/main.nf',
          '-c', meson.current_build_dir() / 'oncoscanner' / 'user.config',
          '--workflow_variation', 'clc',
          '--samplesheet', test_datadir/'test_clc_dna_only_samplesheet.csv',
          '--clc_workflow_name' ,clc_workflow_name,
          '--clc_import_dir' ,clc_import_dir,
          '--clc_export_dir' ,clc_export_dir,
          '--clc_destdir' ,clc_destdir,
          '-with-conda',
          ],
          depends: [clc_nextflow, oncoscanner, create_testdirs, copy_sources, copy_nextpipe_sources, sequence_alignment_nextflow, variantinterpretation, biomarker],
          timeout: 172800 # 2 days
          )

  test('nextflow_workflow_clc_dna_rna_full',
          nextflow_bin,
          is_parallel : false,
          suite : ['nextflow', 'workflow', 'full'],
          env : env_vars,
          args : [
          'run',
          meson.current_build_dir()+'/oncoscanner/main.nf',
          '-c', meson.current_build_dir() / 'oncoscanner' / 'user.config',
          '--workflow_variation', 'clc',
          '--samplesheet', test_datadir/'test_clc_dna_rna_samplesheet.csv',
          '--has_rna true',
          '-with-conda',
          ],
          depends: [clc_nextflow, oncoscanner, create_testdirs, copy_sources, copy_nextpipe_sources, sequence_alignment_nextflow, variantinterpretation],
          timeout: 172800 # 2 days
          )
endif
