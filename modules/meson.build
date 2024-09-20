nextflow_bin = find_program('nextflow')
nextflow_libdir = '/usr/local/lib/nextflow'
test_datadir = get_option('test_datadir')

subdir('arriba_nextflow')
subdir('clc_nextflow')
subdir('sequence_alignment')
subdir('ukb_main_workflow')

nextpipe_sources = [
'nextpipe/resources/usr/bin/organize_variant_table.py',
'nextpipe/resources/usr/bin/variantlist_utils.py',
]

sources = ['ukb_main_workflow'/'test_main.config']
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
	output : ['organize_variant_table.py','variantlist_utils.py'],
	command : ['mkdir', '-p', '@OUTDIR@/nextpipe/resources/usr/bin', '&&', 'cp', '-r', '@INPUT0@', '@INPUT1@', '@OUTDIR@/nextpipe/resources/usr/bin/'],
	install : false,
	depends : [create_nextpipe_dir]
	)

copy_nextpipe_sources = custom_target('copy_nextpipe_main',
	input : 'nextpipe/main.nf',
	output : 'main.nf',
	command : ['cp', '-r', meson.current_source_dir() / '@INPUT@', '@OUTDIR@/nextpipe/'],
	install : false,
	depends : [copy_nextpipe_bin]
	)


copy_main_sources = custom_target('copy_main_sources',
	input : sources + testdata,
	output : 'test_main.config',
	command : ['cp', '-r', '@INPUT@', '@OUTDIR@/ukb_main_workflow'],
	install : false
	)

variantinterpretation = custom_target('variantinterpretation',
	input : 'variantinterpretation/main.nf',
	output : 'variantinterpretation',
	command : ['cp', '-r', '@CURRENT_SOURCE_DIR@/variantinterpretation', '@OUTPUT@/'],
	depends: [copy_main_sources, arriba_nextflow, clc_nextflow]
	)



test('nextflow_workflow_variantinterpretation',
	nextflow_bin,
	is_parallel : false,
	suite : ['nextflow', 'workflow'],
	env: ['NEXTFLOW_MODULES='+meson.current_build_dir(), 
		'NAS_IMPORT_DIR='+meson.current_build_dir()/'nas_import', 
		'NAS_EXPORT_DIR='+meson.current_build_dir()/'nas_export'],
	args : [
	'run',
	meson.current_build_dir()+'/ukb_main_workflow/main.nf',
	'-c', meson.current_build_dir()/ 'ukb_main_workflow' / 'test_main.config',
	'--workflow_variation', 'variantinterpretation',
	'--samplesheet', test_datadir / 'test_variantinterpretation_samplesheet.csv',
	],
	depends: [copy_main_sources, variantinterpretation, ukb_main_workflow],
	timeout: 6000
	)

test('nextflow_workflow_arriba',
	nextflow_bin,
	is_parallel : false,
	suite : ['nextflow', 'workflow'],
	env: ['NEXTFLOW_MODULES='+meson.current_build_dir(), 
		'NAS_IMPORT_DIR='+meson.current_build_dir()/'nas_import', 
		'NAS_EXPORT_DIR='+meson.current_build_dir()/'nas_export'],
	args : [
	'run',
	meson.current_build_dir()+'/ukb_main_workflow/main.nf',
	'-c', meson.current_build_dir()/ 'ukb_main_workflow' / 'test_main.config',
	'--workflow_variation', 'arriba',
	'--samplesheet', test_datadir / 'test_arriba_samplesheet.csv',
	],
	depends: [copy_main_sources, arriba_nextflow, ukb_main_workflow],
	timeout: 600
	)


test('nextflow_workflow_sequence_alignment',
	nextflow_bin,
	is_parallel : false,
	suite : ['nextflow', 'workflow'],
	env: ['NEXTFLOW_MODULES='+meson.current_build_dir(), 
		'NAS_IMPORT_DIR='+meson.current_build_dir()/'nas_import', 
		'NAS_EXPORT_DIR='+meson.current_build_dir()/'nas_export'],
	args : [
	'run',
	meson.current_build_dir()+'/ukb_main_workflow/main.nf',
	'-c', meson.current_build_dir()/'ukb_main_workflow' / 'test_main.config',
	'-profile', 'standard',
	'--workflow_variation', 'sequence_alignment',
	'--samplesheet', test_datadir/'test_sequence_alignment_samplesheet.csv',
	],
	depends: [copy_main_sources, variantinterpretation, sequence_alignment_nextflow, ukb_main_workflow],
	timeout: 360000
	)

test('nextflow_workflow_sequence_alignment_local',
	nextflow_bin,
	is_parallel : false,
	suite : ['nextflow', 'workflow'],
	env: ['NEXTFLOW_MODULES='+meson.current_build_dir(), 
		'NAS_IMPORT_DIR='+meson.current_build_dir()/'nas_import', 
		'NAS_EXPORT_DIR='+meson.current_build_dir()/'nas_export'],
	args : [
	'run',
	meson.current_build_dir()+'/ukb_main_workflow/main.nf',
	'-c', meson.current_build_dir()/'ukb_main_workflow' / 'test_main.config',
	'-profile', 'local_nas_workdir',
	'--workflow_variation', 'sequence_alignment',
	'--samplesheet', test_datadir/'test_sequence_alignment_samplesheet.csv',
	],
	depends: [copy_main_sources, variantinterpretation, sequence_alignment_nextflow, ukb_main_workflow],
	timeout: 360000
	)

test('nextflow_workflow_sequence_alignment_hpc',
	nextflow_bin,
	is_parallel : false,
	suite : ['nextflow', 'workflow'],
	env: ['NEXTFLOW_MODULES='+meson.current_build_dir(), 
		'NAS_IMPORT_DIR='+meson.current_build_dir()/'nas_import', 
		'NAS_EXPORT_DIR='+meson.current_build_dir()/'nas_export'],
	args : [
	'run',
	meson.current_build_dir()+'/ukb_main_workflow/main.nf',
	'-c', meson.current_build_dir()/'ukb_main_workflow' / 'test_main.config',
	'-profile', 'hpc',
	'--workflow_variation', 'sequence_alignment',
	'--samplesheet', test_datadir/'test_sequence_alignment_samplesheet.csv',
	],
	depends: [copy_main_sources, variantinterpretation, sequence_alignment_nextflow, ukb_main_workflow],
	timeout: 360000
	)

test('nextflow_workflow_sequence_align_interpret',
	nextflow_bin,
	is_parallel : false,
	suite : ['nextflow', 'workflow'],
	env: ['NEXTFLOW_MODULES='+meson.current_build_dir(), 
		'NAS_IMPORT_DIR='+meson.current_build_dir()/'nas_import', 
		'NAS_EXPORT_DIR='+meson.current_build_dir()/'nas_export'],
	args : [
	'run',
	meson.current_build_dir()+'/ukb_main_workflow/main.nf',
	'-c', meson.current_build_dir()/'ukb_main_workflow' / 'test_main.config',
	'-profile', 'standard',
	'--workflow_variation', 'align_interpret',
	'--samplesheet', test_datadir/'test_sequence_alignment_samplesheet.csv',
	],
	depends: [copy_main_sources, variantinterpretation, sequence_alignment_nextflow, ukb_main_workflow],
	timeout: 360000
	)


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
	env: ['NEXTFLOW_MODULES='+meson.current_build_dir(), 
		'NAS_IMPORT_DIR='+meson.current_build_dir()/'nas_import', 
		'NAS_EXPORT_DIR='+meson.current_build_dir()/'nas_export'],
	args : [
	'run',
	meson.current_build_dir()+'/ukb_main_workflow/main.nf',
	'-c', meson.current_build_dir() / 'ukb_main_workflow' / 'test_main.config',
	'--workflow_variation', 'clc',
	'--samplesheet', test_datadir/'test_clc_samplesheet.csv',
	'-with-conda',
	'-stub'
	],
	depends: [clc_nextflow, ukb_main_workflow, create_testdirs, copy_main_sources, copy_nextpipe_sources, sequence_alignment_nextflow, variantinterpretation],
	timeout: 6000
	)

test('nextflow_workflow_clc_dna_only_full',
	nextflow_bin,
	is_parallel : false,
	suite : ['nextflow', 'workflow', 'full'],
	env: ['NEXTFLOW_MODULES='+meson.current_build_dir(), 
		'NAS_IMPORT_DIR='+meson.current_build_dir()/'nas_import', 
		'NAS_EXPORT_DIR='+meson.current_build_dir()/'nas_export'],
	args : [
	'run',
	meson.current_build_dir()+'/ukb_main_workflow/main.nf',
	'-c', test_datadir/ 'test_main.config',
	'--workflow_variation', 'clc',
	'--samplesheet', test_datadir/'test_clc_dna_only_samplesheet.csv',
	'-with-conda',
	],
	depends: [clc_nextflow, ukb_main_workflow, create_testdirs, copy_main_sources, copy_nextpipe_sources, sequence_alignment_nextflow, variantinterpretation],
	timeout: 172800 # 2 days
	)

test('nextflow_workflow_clc_dna_rna_full',
	nextflow_bin,
	is_parallel : false,
	suite : ['nextflow', 'workflow', 'full'],
	env: ['NEXTFLOW_MODULES='+meson.current_build_dir(), 
		'NAS_IMPORT_DIR='+meson.current_build_dir()/'nas_import', 
		'NAS_EXPORT_DIR='+meson.current_build_dir()/'nas_export'],
	args : [
	'run',
	meson.current_build_dir()+'/ukb_main_workflow/main.nf',
	'-c', test_datadir/ 'test_main.config',
	'--workflow_variation', 'clc',
	'--samplesheet', test_datadir/'test_clc_dna_rna_samplesheet.csv',
	'--has_rna true',
	'-with-conda',
	],
	depends: [clc_nextflow, ukb_main_workflow, create_testdirs, copy_main_sources, copy_nextpipe_sources, sequence_alignment_nextflow, variantinterpretation],
	timeout: 172800 # 2 days
	)
