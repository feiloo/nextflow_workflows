nextflow_bin = find_program('nextflow')
nextflow_libdir = '/usr/local/lib/nextflow'
test_datadir = get_option('test_datadir')

subdir('arriba_nextflow')
subdir('clc_nextflow')
subdir('sequence_alignment')
subdir('ukb_main_workflow')

sources = ['ukb_main_workflow'/'test_main.config']
testdata = []

copy_sources = custom_target('copy_sources',
	input : sources + testdata,
	output : 'test_main.config',
	command : ['cp', '-r', '@INPUT@', '@OUTDIR@/ukb_main_workflow'],
	install : false
	)

variantinterpretation = custom_target('variantinterpretation',
	input : 'variantinterpretation/main.nf',
	output : 'variantinterpretation',
	command : ['cp', '-r', '@CURRENT_SOURCE_DIR@/variantinterpretation', '@OUTPUT@/'],
	depends: [copy_sources, arriba_nextflow, clc_nextflow]
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
	depends: [copy_sources, variantinterpretation, ukb_main_workflow],
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
	depends: [copy_sources, arriba_nextflow, ukb_main_workflow],
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
	depends: [copy_sources, variantinterpretation, sequence_alignment_nextflow, ukb_main_workflow],
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
	depends: [copy_sources, variantinterpretation, sequence_alignment_nextflow, ukb_main_workflow],
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
	depends: [copy_sources, variantinterpretation, sequence_alignment_nextflow, ukb_main_workflow],
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
	depends: [copy_sources, variantinterpretation, sequence_alignment_nextflow, ukb_main_workflow],
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
	'-stub'
	],
	depends: [clc_nextflow, ukb_main_workflow, create_testdirs, copy_sources, sequence_alignment_nextflow, variantinterpretation],
	timeout: 6000
	)
