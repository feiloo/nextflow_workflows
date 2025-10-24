install_subdir('../oncoscanner', 
	install_dir:nextflow_libdir, 
	exclude_files:['meson.build'], 
	install_tag : 'nextflow')

nextflow_files = ['main.nf','nextflow.config','user.config']


if enable_proprietary
  deps = [arriba_nextflow, clc_nextflow, sequence_alignment_nextflow, variantinterpretation]
else
  deps = [arriba_nextflow, sequence_alignment_nextflow, variantinterpretation]
endif

oncoscanner = custom_target('oncoscanner',
	input : nextflow_files,
	output : '.',
	command : ['cp', '@INPUT@', '@OUTPUT@/'],
	depends: deps
	)
