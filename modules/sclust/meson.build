install_subdir('../sclust', 
	install_dir:nextflow_libdir, 
	exclude_files:['meson.build'], 
	install_tag : 'nextflow')

sclust = custom_target('sclust',
	input : 'main.nf',
	output : 'main.nf',
	command : ['cp', '@INPUT@', '@OUTPUT@']
	)
