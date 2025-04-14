install_subdir('../biomarker', 
	install_dir:nextflow_libdir, 
	exclude_files:['meson.build'], 
	install_tag : 'nextflow')

nextflow_files = ['main.nf']

biomarker = custom_target('biomarker',
	input : nextflow_files,
	output : '.',
	command : ['cp', '@INPUT@', '@OUTPUT@/'],
	)
