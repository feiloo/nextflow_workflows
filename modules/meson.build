nextflow_libdir = '/usr/local/lib/nextflow'
testconfig_path = '/home/fhoelsch/nextflow_conf_specific.config'

subdir('arriba_nextflow')
subdir('clc_nextflow')


variantinterpretation = custom_target('variantinterpretation',
	input : 'variantinterpretation/main.nf',
	output : 'variantinterpretation',
	command : ['cp', '-r', '@CURRENT_SOURCE_DIR@/variantinterpretation', '@OUTPUT@/'],
	depends: [arriba_nextflow, clc_nextflow]
	)

subdir('ukb_main_workflow')

nextflow_bin = find_program('nextflow')

test('nextflow_workflow',
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
	timeout: 600
	)
