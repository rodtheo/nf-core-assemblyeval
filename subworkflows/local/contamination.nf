include { FCS_FCSADAPTOR } from '../../modules/nf-core/fcs/fcsadaptor/main'
include { FCS_FCSGX } from '../../modules/nf-core/fcs/fcsgx/main'
include { FCS_CLEAN } from '../../modules/local/fcs/fcsclean/main'


workflow CONTAMINATION_ASM {

	take:
    reads        // queue channel; [ sample_id, [ file(read1), file(read2) ] ]
    reference    // file( "path/to/reference" )

    main:
    // Quality Check input reads
    // READ_QC ( reads )

    // Align reads to reference
    Channel.empty()
        .set { results_dir_ch }

	Channel.empty()
        .set { results_summary_ch }

	Channel.empty()
		.set { version_ch }

	Channel.empty()
		.set { cleaned_fasta_ch }

	reference.map{ meta, assemblies -> {
		def meta_new = meta.clone();
		meta_new.put("id", meta.get("id")+"_"+meta.get("build"));
		meta_new.put("reference_id", meta.get("id"));
		[ meta_new, assemblies['pri_asm'] ]
	 } }.set{ reference_ch }

	FCS_FCSADAPTOR( reference_ch )

	version_ch.mix(FCS_FCSADAPTOR.out.versions).set{ version_ch }

	// dbgx_ch = Channel.fromPath('./data/fcs-database')
	dbgx_ch = Channel.fromPath('/home/rodtheo/Bioinfo/data/fcs-database')

	reference_ch.combine( dbgx_ch ).multiMap{ meta, asm, db_path ->
		assembly: [meta, asm]
		db: [db_path] }.set{ to_fcsgx_ch }

	FCS_FCSGX( to_fcsgx_ch.assembly, to_fcsgx_ch.db )

	version_ch.mix(FCS_FCSGX.out.versions).set{ version_ch }

	reference_ch.combine(FCS_FCSGX.out.fcs_gx_report, by: [0]).set{ to_clean_fa_ch }

	to_clean_fa_ch.view{ "\n########### TO CLEAN: $it\n" }

	FCS_CLEAN( to_clean_fa_ch )

	// cleaned_fasta_ch.mix( FCS_CLEAN.out.fcs_clean_fasta )

	FCS_CLEAN.out.fcs_clean_fasta.map{ meta, asm_cleaned -> 
		def meta_new = meta.clone();
		meta_new.put("id", meta.get("reference_id"));
		meta_new.put("build", meta.get("build")+"_cleaned");
		[ meta_new, ['pri_asm': asm_cleaned ] ]
	 }.ifEmpty([]).set{ cleaned_fasta_ch }

	out_reference = reference.mix(cleaned_fasta_ch).filter{ it != [] }

	cleaned_fasta_ch.filter{ it != [] }.view{ "CLEANED FASTA: $it" }	

    emit:
	fcs = FCS_FCSGX.out.fcs_gx_report
	asms = out_reference
	versions = version_ch
}