nextflow_libdir = '/usr/local/lib/nextflow'
test_configs = get_option('test_configs')
test_samplesheets = get_option('test_samplesheets')

subdir('arriba_nextflow')
subdir('clc_nextflow')
subdir('sequence_alignment')


variantinterpretation = custom_target('variantinterpretation',
	input : 'variantinterpretation/main.nf',
	output : 'variantinterpretation',
	command : ['cp', '-r', '@CURRENT_SOURCE_DIR@/variantinterpretation', '@OUTPUT@/'],
	depends: [arriba_nextflow, clc_nextflow]
	)

subdir('ukb_main_workflow')

nextflow_bin = find_program('nextflow')

test('nextflow_workflow_variantinterpretation',
	nextflow_bin,
	is_parallel : false,
	suite : ['nextflow', 'workflow'],
	env: ['NEXTFLOW_MODULES='+meson.current_build_dir()],
	args : [
	'run',
	meson.current_build_dir()+'/ukb_main_workflow/main.nf',
	'-c', test_configs/'test_main.config',
	'--workflow_variation', 'variantinterpretation',
	'--samplesheet', test_samplesheets + '/variantinterpretation.csv',
	],
	depends: [variantinterpretation, ukb_main_workflow],
	timeout: 6000
	)

test('nextflow_workflow_arriba',
	nextflow_bin,
	is_parallel : false,
	suite : ['nextflow', 'workflow'],
	env: ['NEXTFLOW_MODULES='+meson.current_build_dir()],
	args : [
	'run',
	meson.current_build_dir()+'/ukb_main_workflow/main.nf',
	'-c', test_configs/'test_main.config',
	'--workflow_variation', 'arriba',
	],
	depends: [arriba_nextflow, ukb_main_workflow],
	timeout: 600
	)

test('nextflow_workflow_sequence_alignment',
	nextflow_bin,
	is_parallel : false,
	suite : ['nextflow', 'workflow'],
	env: ['NEXTFLOW_MODULES='+meson.current_build_dir()],
	args : [
	'run',
	meson.current_build_dir()+'/ukb_main_workflow/main.nf',
	'-c', test_configs/'test_main.config',
	'--workflow_variation', 'sequence_alignment',
	'--samplesheet', test_samplesheets/'test_sequence_alignment_samplesheet.csv',
	],
	depends: [variantinterpretation, sequence_alignment_nextflow, ukb_main_workflow],
	timeout: 3600
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
	'-c', test_configs/'test_main.config',
	'--workflow_variation', 'clc',
	'--samplesheet', test_samplesheets/'test_clc_samplesheet.csv',
	'-stub'
	],
	depends: [clc_nextflow, ukb_main_workflow, create_testdirs],
	timeout: 6000
	)

