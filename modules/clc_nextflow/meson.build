container_engine = find_program('podman')

install_subdir('../clc_nextflow', 
	install_dir:nextflow_libdir, 
	exclude_files:['meson.build'], 
	install_tag : 'nextflow')

clc_nextflow = custom_target('clc_nextflow',
	input : 'main.nf',
	output : 'main.nf',
	command : ['cp', '@INPUT@', '@OUTPUT@']
	)


clc_cli_tools_file = 'CLCServerCommandLineTools_24_0_1_64.sh'
containerfile = 'clc_client.containerfile'
container_image = 'clc_client.tar'

sources = ['third_party'/clc_cli_tools_file]
testdata = []

copy_sources = custom_target('copy_sources',
	input : sources + testdata,
	output : clc_cli_tools_file,
	command : ['cp', '-r', '@INPUT@', '@OUTDIR@'],
	install : false
	)

clc_nextflow_container = custom_target('clc_nextflow_container',
	input : [containerfile, 'third_party'/clc_cli_tools_file],
	output : container_image,
	command : [container_engine, 'build', '-t', 
		'clc_client:latest', 
		'--build-arg', 'clc_servertools_install_sh='+clc_cli_tools_file,
		'-f', '@INPUT0@', '-o', 'type=tar,dest=@OUTPUT@', meson.current_build_dir()],
	build_by_default : false,
	depends : [copy_sources],
	install : false
	)

