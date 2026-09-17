install_subdir('../fusions', 
	install_dir:nextflow_libdir, 
	exclude_files:['meson.build'], 
	install_tag : 'nextflow')

arriba_nextflow = custom_target('fusions',
	input : 'main.nf',
	output : 'main.nf',
	command : ['cp', '@INPUT@', '@OUTPUT@']
	)
