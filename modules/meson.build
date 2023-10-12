nextflow_libdir = '/usr/local/lib/nextflow'
testconfig_path = '/home/fhoelsch/nextflow_conf_specific.conf'

subdir('arriba_nextflow')
subdir('clc_nextflow')
subdir('ukb_main_workflow')

nextflow_bin = find_program('nextflow')

test('nextflow_workflow',
	nextflow_bin,
	is_parallel : false,
	suite : ['nextflow', 'workflow'],
	args : [
	'run',
	meson.current_build_dir()+'ukb_main_workflow/main.nf',
	'-c', testconfig_path,
	],
	timeout: 600
	)
