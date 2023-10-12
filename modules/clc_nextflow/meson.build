install_subdir('../clc_nextflow', 
	install_dir:nextflow_libdir, 
	exclude_files:['meson.build'], 
	install_tag : 'nextflow')

custom_target('clc_nextflow',
	input : 'main.nf',
	output : 'main.nf',
	command : ['cp', '@INPUT@', '@OUTPUT@']
	)
