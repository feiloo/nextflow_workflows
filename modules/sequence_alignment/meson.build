install_subdir('../sequence_alignment', 
	install_dir:nextflow_libdir, 
	exclude_files:['meson.build'], 
	install_tag : 'nextflow')

nextflow_files = ['main.nf','nextflow.config', 'fastqc.nf', 'fastp.nf', 'bwamem2.nf',
	'align_and_postprocess.nf','bwa.nf','gatk_mutect.nf',
	'gatk.nf', 'samtools.nf', 'utils.nf'
	]
sequence_alignment_nextflow = custom_target('sequence_alignment_nextflow',
	input : nextflow_files,
	output : '.',
	command : ['cp', '@INPUT@', '@OUTPUT@/'],
	depends: []
	)
