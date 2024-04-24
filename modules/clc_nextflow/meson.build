container_engine = find_program('podman')
external_deps = get_option('external_deps')

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

sources = ['download_clc_cli_tools.sh']
testdata = []

copy_sources = custom_target('copy_sources',
	input : sources + testdata,
	output : sources[0],
	command : ['cp', '-r', '@INPUT@', '@OUTDIR@'],
	install : false
	)

download_clc_tools = custom_target('download_clc_tools',
	input : 'download_clc_cli_tools.sh',
	output : clc_cli_tools_file,
	command : ['download_clc_cli_tools.sh', external_deps, clc_cli_tools_file, '&&', 'cp', external_deps / clc_cli_tools_file, meson.current_build_dir() / clc_cli_tools_file],
	depends: [],
	install : false,
	build_by_default : false
	)

clc_nextflow_container = custom_target('clc_nextflow_container',
	input : [containerfile],
	output : container_image,
	command : [container_engine, 'build', '-t', 
		'clc_client:latest', 
		'--build-arg', 'clc_servertools_install_sh='+clc_cli_tools_file,
		'-f', '@INPUT0@', '-o', 'type=tar,dest=@OUTPUT@', meson.current_build_dir()],
	build_by_default : false,
	depends : [copy_sources],
	install : false
	)

