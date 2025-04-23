nextflow_bin = find_program('nextflow')
nextflow_libdir = '/usr/local/lib/nextflow'
test_datadir = get_option('test_datadir')
enable_proprietary = get_option('enable_proprietary')

subdir('arriba_nextflow')

nextpipe_sources = [
'nextpipe/resources/usr/bin/transcript_list_08_08_2024.py',
'nextpipe/resources/usr/bin/organize_variant_table.py',
'nextpipe/resources/usr/bin/variantlist_utils.py',
]

if enable_proprietary
  subdir('clc_nextflow')
endif

subdir('sequence_alignment')
subdir('biomarker')

sources = ['ukb_main_workflow'/'user.config']
testdata = []

create_nextpipe_dir = custom_target('create_nextpipe_dir',
	command : ['mkdir', '-p', meson.current_build_dir() / 'nextpipe'/'bin', 
		'&&', 
		'touch', 'nextpipe_testfile'
		],
	output : ['nextpipe_testfile']
	)

copy_nextpipe_bin = custom_target('copy_nextpipe_bin',
	input : nextpipe_sources,
	output : ['organize_variant_table.py','variantlist_utils.py', 'transcript_list_08_08_2024'],
	command : ['mkdir', '-p', '@OUTDIR@/nextpipe/resources/usr/bin', '&&', 'cp', '-r', '@INPUT@', '@OUTDIR@/nextpipe/resources/usr/bin/'],
	install : false,
	depends : [create_nextpipe_dir]
	)

copy_nextpipe_sources = custom_target('copy_nextpipe_main',
	input : 'nextpipe/main.nf',
	output : 'main.nf',
	command : ['cp', '-r', '@INPUT@', '@OUTDIR@/nextpipe/'],
	install : false,
	depends : [copy_nextpipe_bin]
	)


copy_main_sources = custom_target('copy_main_sources',
	input : sources + testdata,
	output : 'user.config',
	command : ['cp', '-r', '@INPUT@', '@OUTDIR@/ukb_main_workflow'],
	install : false
	)

# collected by running "find" for .nf environment.yml meta.yml .sh .py and .config files
variantinterpretation_files = [
  'variantinterpretation/modules/local/datavzrd/datavzrd/main.nf',
  'variantinterpretation/modules/local/datavzrd/preprocess/main.nf',
  'variantinterpretation/modules/local/gatk4/validatevariants/main.nf',
  'variantinterpretation/modules/local/multiqcreport_warnings.nf',
  'variantinterpretation/modules/local/tmbcalculation/main.nf',
  'variantinterpretation/modules/local/vcfchecks/main.nf',
  'variantinterpretation/modules/local/vembrane/create-vembrane-fields/main.nf',
  'variantinterpretation/modules/local/vembrane/create-vembrane-tags/main.nf',
  'variantinterpretation/modules/local/vembrane/filter/main.nf',
  'variantinterpretation/modules/local/vembrane/table/main.nf',
  'variantinterpretation/modules/local/vembrane/tag/main.nf',
  'variantinterpretation/modules/local/bcftools/samplerename/main.nf',
  'variantinterpretation/modules/local/bedfile/checkbedfile/main.nf',
  'variantinterpretation/modules/local/bedfile/preparebedfile/main.nf',
  'variantinterpretation/modules/nf-core/multiqc/main.nf',
  'variantinterpretation/modules/nf-core/bcftools/index/main.nf',
  'variantinterpretation/modules/nf-core/bcftools/norm/main.nf',
  'variantinterpretation/modules/nf-core/bcftools/stats/main.nf',
  'variantinterpretation/modules/nf-core/bcftools/view/main.nf',
  'variantinterpretation/modules/nf-core/bcftools/annotate/main.nf',
  'variantinterpretation/modules/nf-core/bcftools/merge/main.nf',
  'variantinterpretation/modules/nf-core/samtools/dict/main.nf',
  'variantinterpretation/modules/nf-core/samtools/faidx/main.nf',
  'variantinterpretation/modules/nf-core/ensemblvep/filtervep/main.nf',
  'variantinterpretation/modules/nf-core/ensemblvep/vep/main.nf',
  'variantinterpretation/modules/nf-core/tabix/bgzip/main.nf',
  'variantinterpretation/modules/nf-core/tabix/tabix/main.nf',
  'variantinterpretation/subworkflows/local/html_report/main.nf',
  'variantinterpretation/subworkflows/local/variantfilter/main.nf',
  'variantinterpretation/subworkflows/local/check_vcf/main.nf',
  'variantinterpretation/subworkflows/local/merge_vcfs/main.nf',
  'variantinterpretation/subworkflows/local/process_vcf/main.nf',
  'variantinterpretation/subworkflows/local/tsv_conversion/main.nf',
  'variantinterpretation/subworkflows/local/utils_nfcore_variantinterpretation_pipeline/main.nf',
  'variantinterpretation/subworkflows/local/vcf_roi_tagging/main.nf',
  'variantinterpretation/subworkflows/nf-core/utils_nextflow_pipeline/main.nf',
  'variantinterpretation/subworkflows/nf-core/utils_nfcore_pipeline/main.nf',
  'variantinterpretation/subworkflows/nf-core/utils_nfvalidation_plugin/main.nf',
  'variantinterpretation/workflows/variantinterpretation.nf',
  'variantinterpretation/main.nf',
  'variantinterpretation/modules/local/bcftools/samplerename/environment.yml',
  'variantinterpretation/modules/nf-core/multiqc/environment.yml',
  'variantinterpretation/modules/nf-core/bcftools/index/environment.yml',
  'variantinterpretation/modules/nf-core/bcftools/norm/environment.yml',
  'variantinterpretation/modules/nf-core/bcftools/stats/environment.yml',
  'variantinterpretation/modules/nf-core/bcftools/view/environment.yml',
  'variantinterpretation/modules/nf-core/bcftools/annotate/environment.yml',
  'variantinterpretation/modules/nf-core/bcftools/merge/environment.yml',
  'variantinterpretation/modules/nf-core/samtools/dict/environment.yml',
  'variantinterpretation/modules/nf-core/samtools/faidx/environment.yml',
  'variantinterpretation/modules/nf-core/ensemblvep/filtervep/environment.yml',
  'variantinterpretation/modules/nf-core/ensemblvep/vep/environment.yml',
  'variantinterpretation/modules/nf-core/tabix/bgzip/environment.yml',
  'variantinterpretation/modules/nf-core/tabix/tabix/environment.yml',
  'variantinterpretation/modules/local/datavzrd/preprocess/resources/usr/bin/preprocess_datavzrd.py',
  'variantinterpretation/modules/local/tmbcalculation/resources/usr/bin/calculate_TMB.py',
  'variantinterpretation/modules/local/vcfchecks/resources/usr/bin/check_vcf.py',
  'variantinterpretation/modules/local/vembrane/create-vembrane-fields/resources/usr/bin/create_vembrane_fields.py',
  'variantinterpretation/modules/local/vembrane/create-vembrane-tags/resources/usr/bin/create_vembrane_tags.py',
  'variantinterpretation/modules/local/resources/bin/check_bedfiles.py',
  'variantinterpretation/modules/nf-core/custom/dumpsoftwareversions/templates/dumpsoftwareversions.py',
  'variantinterpretation/bin/check_vcf.py',
  'variantinterpretation/bin/create_vembrane_fields.py',
  'variantinterpretation/bin/calculate_TMB.py',
  'variantinterpretation/bin/check_bedfiles.py',
  'variantinterpretation/bin/create_vembrane_tags.py',
  'variantinterpretation/bin/preprocess_datavzrd.py',
  ]

# just builds/copies the variantinterpretation, since it doesnt have meson
variantinterpretation = custom_target('variantinterpretation',
	input : variantinterpretation_files,
	output : 'variantinterpretation',
	command : ['cp', '-rT', '@CURRENT_SOURCE_DIR@/variantinterpretation', '@OUTPUT@/'],
	)

# need to include here, so that variantinterpretation is defined
subdir('ukb_main_workflow')

env_vars = ['NEXTFLOW_MODULES='+meson.current_build_dir(), 
              'NAS_IMPORT_DIR='+meson.current_build_dir()/'nas_import', 
              'NAS_EXPORT_DIR='+meson.current_build_dir()/'nas_export']

test('nextflow_workflow_variantinterpretation',
	nextflow_bin,
	is_parallel : false,
	suite : ['nextflow', 'workflow'],
	env: env_vars,
	args : [
	'run',
	meson.current_build_dir()+'/ukb_main_workflow/main.nf',
	'-c', meson.current_build_dir()/ 'ukb_main_workflow' / 'user.config',
	'--workflow_variation', 'variantinterpretation',
	'--samplesheet', test_datadir / 'test_variantinterpretation_samplesheet.csv',
	],
	depends: [copy_main_sources, variantinterpretation, ukb_main_workflow],
	timeout: 6000
	)

test('nextflow_workflow_arriba',
	nextflow_bin,
	is_parallel : false,
	suite : ['nextflow', 'workflow'],
	env: env_vars,
	args : [
	'run',
	meson.current_build_dir()+'/ukb_main_workflow/main.nf',
	'-c', meson.current_build_dir()/ 'ukb_main_workflow' / 'user.config',
	'--workflow_variation', 'arriba',
	'--samplesheet', test_datadir / 'test_arriba_samplesheet.csv',
	],
	depends: [copy_main_sources, arriba_nextflow, ukb_main_workflow],
	timeout: 600
	)


test('nextflow_workflow_sequence_alignment',
	nextflow_bin,
	is_parallel : false,
	suite : ['nextflow', 'workflow'],
	env: env_vars,
	args : [
	'run',
	meson.current_build_dir()+'/ukb_main_workflow/main.nf',
	'-c', meson.current_build_dir()/'ukb_main_workflow' / 'user.config',
	'-profile', 'standard',
	'--workflow_variation', 'sequence_alignment',
	'--samplesheet', test_datadir/'test_sequence_alignment_samplesheet.csv',
	],
	depends: [copy_main_sources, variantinterpretation, sequence_alignment_nextflow, ukb_main_workflow],
	timeout: 360000
	)

test('nextflow_workflow_sequence_alignment_local',
	nextflow_bin,
	is_parallel : false,
	suite : ['nextflow', 'workflow'],
	env: env_vars,
	args : [
	'run',
	meson.current_build_dir()+'/ukb_main_workflow/main.nf',
	'-c', meson.current_build_dir()/'ukb_main_workflow' / 'user.config',
	'-profile', 'local_nas_workdir',
	'--workflow_variation', 'sequence_alignment',
	'--samplesheet', test_datadir/'test_sequence_alignment_samplesheet.csv',
	],
	depends: [copy_main_sources, variantinterpretation, sequence_alignment_nextflow, ukb_main_workflow],
	timeout: 360000
	)

test('nextflow_workflow_sequence_alignment_hpc',
	nextflow_bin,
	is_parallel : false,
	suite : ['nextflow', 'workflow'],
	env: env_vars,
	args : [
	'run',
	meson.current_build_dir()+'/ukb_main_workflow/main.nf',
	'-c', meson.current_build_dir()/'ukb_main_workflow' / 'user.config',
	'-profile', 'hpc',
	'--workflow_variation', 'sequence_alignment',
	'--samplesheet', test_datadir/'test_sequence_alignment_samplesheet.csv',
	],
	depends: [copy_main_sources, variantinterpretation, sequence_alignment_nextflow, ukb_main_workflow],
	timeout: 360000
	)

test('nextflow_workflow_sequence_align_interpret',
	nextflow_bin,
	is_parallel : false,
	suite : ['nextflow', 'workflow'],
	env: env_vars,
	args : [
	'run',
	meson.current_build_dir()+'/ukb_main_workflow/main.nf',
	'-c', meson.current_build_dir()/'ukb_main_workflow' / 'user.config',
	'-profile', 'standard',
	'--workflow_variation', 'align_interpret',
	'--samplesheet', test_datadir/'test_sequence_alignment_samplesheet.csv',
	],
	depends: [copy_main_sources, variantinterpretation, sequence_alignment_nextflow, biomarker, ukb_main_workflow],
	timeout: 360000
	)



if enable_proprietary
  create_testdirs = custom_target('create_testdirs',
          command : ['mkdir', '-p', meson.current_build_dir() / 'nas_import', 
                  meson.current_build_dir() / 'nas_export', 
                  '&&', 'touch', 'nas_testfile'
                  ],
          output : ['nas_testfile']
          )

  test('nextflow_workflow_clc',
          nextflow_bin,
          is_parallel : false,
          suite : ['nextflow', 'workflow'],
          env : env_vars,
          args : [
          'run',
          meson.current_build_dir()+'/ukb_main_workflow/main.nf',
          '-c', meson.current_build_dir() / 'ukb_main_workflow' / 'user.config',
          '--workflow_variation', 'clc',
          '--samplesheet', test_datadir/'test_clc_samplesheet.csv',
          '-with-conda',
          '-stub'
          ],
          depends: [clc_nextflow, ukb_main_workflow, create_testdirs, copy_sources, copy_nextpipe_sources, sequence_alignment_nextflow, variantinterpretation],
          timeout: 6000
          )

  test('nextflow_workflow_clc_dna_only_full',
          nextflow_bin,
          is_parallel : false,
          suite : ['nextflow', 'workflow', 'full'],
          env : env_vars,
          args : [
          'run',
          meson.current_build_dir()+'/ukb_main_workflow/main.nf',
          '-c', meson.current_build_dir() / 'ukb_main_workflow' / 'user.config',
          '--workflow_variation', 'clc',
          '--samplesheet', test_datadir/'test_clc_dna_only_samplesheet.csv',
          '-with-conda',
          ],
          depends: [clc_nextflow, ukb_main_workflow, create_testdirs, copy_sources, copy_nextpipe_sources, sequence_alignment_nextflow, variantinterpretation],
          timeout: 172800 # 2 days
          )

  test('nextflow_workflow_clc_dna_rna_full',
          nextflow_bin,
          is_parallel : false,
          suite : ['nextflow', 'workflow', 'full'],
          env : env_vars,
          args : [
          'run',
          meson.current_build_dir()+'/ukb_main_workflow/main.nf',
          '-c', meson.current_build_dir() / 'ukb_main_workflow' / 'user.config',
          '--workflow_variation', 'clc',
          '--samplesheet', test_datadir/'test_clc_dna_rna_samplesheet.csv',
          '--has_rna true',
          '-with-conda',
          ],
          depends: [clc_nextflow, ukb_main_workflow, create_testdirs, copy_sources, copy_nextpipe_sources, sequence_alignment_nextflow, variantinterpretation],
          timeout: 172800 # 2 days
          )
endif
