nextflow_libdir = '/usr/local/lib/nextflow'
testconfig_path = get_option('testconfig_path')
test_samplesheet = get_option('test_samplesheet')

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
	'-c', testconfig_path,
	'--workflow_variation', 'variantinterpretation',
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
	'-c', testconfig_path,
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
	'-c', testconfig_path,
	'--workflow_variation', 'sequence_alignment',
	'--samplesheet', test_samplesheet,
	],
	depends: [variantinterpretation, sequence_alignment_nextflow, ukb_main_workflow],
	timeout: 3600
	)
