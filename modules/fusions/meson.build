install_subdir('../fusions', 
	install_dir:nextflow_libdir, 
	exclude_files:['meson.build'], 
	install_tag : 'nextflow')

fusions = custom_target('fusions',
	input : 'main.nf',
	output : 'main.nf',
	command : ['cp', '@INPUT@', '@OUTPUT@']
	)
