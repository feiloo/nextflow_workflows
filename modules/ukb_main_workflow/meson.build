install_subdir('../ukb_main_workflow', 
	install_dir:nextflow_libdir, 
	exclude_files:['meson.build'], 
	install_tag : 'nextflow')

nextflow_files = ['main.nf','nextflow.config']
ukb_main_workflow = custom_target('ukb_main_workflow',
	input : nextflow_files,
	output : '.',
	command : ['cp', '@INPUT@', '@OUTPUT@/'],
	depends: [arriba_nextflow, clc_nextflow]
	)
