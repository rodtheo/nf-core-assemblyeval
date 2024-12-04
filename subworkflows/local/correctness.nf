include { ALE } from '../../modules/nf-core/ale/main'
include { UCSC_WIGTOBIGWIG } from '../../modules/nf-core/ucsc/wigtobigwig/main'
include { ALE_TO_WIGGLE } from '../../modules/local/ALE_towig' 
include { CUSTOM_GETCHROMSIZES } from '../../modules/nf-core/custom/getchromsizes/main'
include { REAPR } from '../../modules/local/REAPR'
include { UCSC_BIGWIGTOBEDGRAPH } from '../../modules/local/bigwigtobedgraph'
include { HEADER_FASTA_REAPR } from '../../modules/local/header_fasta_reapr'
include { IGVREPORTS } from '../../modules/nf-core/igvreports/main'
include { PREPARE_REPORT_IGV } from '../../modules/local/prepare_report_igv'
include { SUBSET_REAPR } from '../../modules/local/subset_reapr'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX; TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_REAPR } from '../../modules/nf-core/tabix/bgziptabix/main'
include { REAPER_BEDGRAPH } from '../../modules/local/REAPR_BEDGRAPH'

workflow CORRECTNESS_ASM { 

    take:
	joined_ch
	asm 
    bam        // queue channel: [ sample_id, file(bam_file) ]
    bai    // queue channel: [ sample_id, file(bai_file) ]

    main:
    // Quality Check input reads
    // READ_QC ( reads )

    
	// reference.map{ meta, assemblies -> {
	// 	def meta_new = meta.clone();
	// 	meta_new.put("id", meta.get("id")+"_"+meta.get("build"));
	// 	[ meta_new, assemblies['pri_asm'] ]
	//  } }.set{ reference_ch }
	

	// bam.combine(bai).map{ meta_bam, bam, meta_bai, bai -> {
	// 	def meta_cor_reads = meta_reads.clone();
	// 	meta_cor_reads.put("id", meta_asm.get("id"));
	// 	[ meta_cor_reads, reads, meta_asm, asm ]
	// }
	// }.set{ inputs_reads_asm_ch }
	// inputs_reads_asm_ch.view{ "ASSEMBLIES+READS COMBINED CHANNEL: $it"}

	// left = asm.map{ meta_asm, asm -> [meta_asm['id'], meta_asm, asm] }
	// right = bam.map{ meta_bam, bam -> [meta_bam['id'], meta_bam, bam] }

	// joined_ch = left.join(right, failOnMismatch: true).view{ "### JOINED: $it"}

	joined_ch.multiMap{ meta_id, asm, bam ->
		asm: [['id': meta_id], asm[1], bam[1]] }.set{ input_ale_ch }

	// SAMTOOLS_SORT( input_samtools_sort_ch.bam,  input_samtools_sort_ch.ref )

	// Align reads to reference
    Channel.empty()
        .set { ale_out_ch }

	ALE ( input_ale_ch )

	ale_out_ch.mix(ALE.out.ale).set{ ale_out_ch }

	Channel.empty()
        .set { reapr_out_ch }

	Channel.empty()
        .set { igv_report_out_ch }

	joined_ch.multiMap{ meta_id, asm, bam ->
		asm: asm
		reapr: bam }.set{ input_reapr_ch }


	HEADER_FASTA_REAPR ( input_reapr_ch.asm )

	HEADER_FASTA_REAPR.out.asm.map{ meta, asm -> [meta['id'], meta, asm]}.set{ left_reapr_ch }
	input_reapr_ch.reapr.map{ meta, bam -> [meta['id'], meta, bam]}.set{ right_reapr_ch }

	left_reapr_ch.join(right_reapr_ch).map{ meta_id, meta_asm, asm, meta_bam, bam -> [meta_asm, asm, bam]}.set{ reapr_ch_header }
	// reapr_ch_header.view{ "\n--------------\nREEEEAPR HEADER: $it"}

	// HEADER_FASTA_REAPR.out.asm.join(input_reapr_ch.reapr, by: [0]).set { reapr_ch_header }

	REAPR ( reapr_ch_header )
	
	// SUB-SUB-WORKFLOW CREATE IGV REPORTS
	joined_ch.multiMap{ meta_id, asm, bam ->
		asm: [['id': meta_id], asm[1]] }.set{ input_chr_size_ch }

	// input_chr_size_ch.view{ "\n--------------\nGET CHR SIZE HEADER: $it \n\n"}

	ALE_TO_WIGGLE( ale_out_ch )


	CUSTOM_GETCHROMSIZES( input_chr_size_ch.asm )

	// ALE_TO_WIGGLE.out.ale_wig_base.view{ "\n--------------\nGET ALE WIGGLE: $it \n\n"}
	// CUSTOM_GETCHROMSIZES.out.sizes.view{ "\n--------------\nGET CHR SIZE HEADER: $it \n\n"}

	ALE_TO_WIGGLE.out.ale_wig_base.set{ left_wiggle }
	CUSTOM_GETCHROMSIZES.out.sizes.set{ right_wiggle }
	ALE_TO_WIGGLE.out.ale_wig_base.join(right_wiggle)
	.join(ALE_TO_WIGGLE.out.ale_wig_depth)
	.join(ALE_TO_WIGGLE.out.ale_wig_insert)
	.join(ALE_TO_WIGGLE.out.ale_wig_kmer)
	.join(ALE_TO_WIGGLE.out.ale_wig_place)
	.multiMap{ meta_id, wig_in, sizes_in, wig_depth, wig_insert, wig_kmer, wig_place  -> 
		wig_base: [[['id': meta_id['id']+'_base','id_original': meta_id['id']], wig_in], [sizes_in]]
		wig_depth: [[['id': meta_id['id']+'_depth','id_original': meta_id['id']], wig_depth], [sizes_in]]
		wig_insert: [[['id': meta_id['id']+'_insert','id_original': meta_id['id']], wig_insert], [sizes_in]]
		wig_kmer: [[['id': meta_id['id']+'_kmer','id_original': meta_id['id']], wig_kmer], [sizes_in]]
		wig_place: [[['id': meta_id['id']+'_place','id_original': meta_id['id']], wig_place], [sizes_in]]
		}.set{ in_wig_ch }

	// in_wig_ch.wig_base.view{ "\n--------------\nWIG TO BEDGRAPH: $it \n\n" }

	in_wig_ch.wig_base.mix(in_wig_ch.wig_depth, in_wig_ch.wig_insert,
	in_wig_ch.wig_kmer, in_wig_ch.wig_place).multiMap{ wig_in, sizes_in -> 
		wig: wig_in
		sizes_in: sizes_in}.set{ wig_base_ch }

	wig_base_ch.wig.view{ "\n--------------\nWIG TO BEDGRAPH: $it \n\n" }

	UCSC_WIGTOBIGWIG( wig_base_ch.wig, wig_base_ch.sizes_in )

	UCSC_BIGWIGTOBEDGRAPH(UCSC_WIGTOBIGWIG.out)

	// def ch_render_igv = true

	SUBSET_REAPR( REAPR.out.reapr_score_errors )

	SUBSET_REAPR.out.reaper_bed_failure.ifEmpty('\n\n&&&&&&&&&&&&&&&&&&&&&&&&\n REAPR EMPTY\n&&&&&&&&&&&&&&&&&&&&&&&\n\n').view()

	SUBSET_REAPR.out.reaper_bed_failure.countLines().branch{ v -> 
		keep_going: v > 0
		no_reaper_failures: v == 0 }.set{ ch_render_igv }

	// if (ch_render_igv) {
	// 	UCSC_WIGTOBIGWIG.out.bw.view{ "\n--------------\nWIG TO BIGWIG: $it \n\n" }
	// }

	Channel.fromPath(params.template_tracks_igv).set{ template_track_igv }

	// template_track_igv.view{ "\n--------------\nTEMPLATE TRACK IGV: $it \n\n"}

	joined_ch.multiMap{ meta_id, asm, bam ->
		meta_id: [meta_id] }.set{ input_chr_meta_id }

	input_chr_meta_id.combine( template_track_igv ).map{ meta_id, template -> [['id': meta_id], template]}.set{ template_track_igv }

	template_track_igv.view{ "\n--------------\nTEMPLATE TRACK UPDATED IGV: $it \n\n"}

	REAPR.out.reapr_dir.map{ meta_complete, reapr_dir -> [['id': meta_complete['id']], reapr_dir]}.set{ reapr_dir_small_meta_ch }

	REAPR.out.reapr_score_per_base.map{ meta_complete, reapr_score_per_base -> [['id': meta_complete['id']], reapr_score_per_base]}.set{ reapr_score_per_base_ch }

	REAPER_BEDGRAPH( reapr_score_per_base_ch )
	TABIX_BGZIPTABIX_REAPR( REAPER_BEDGRAPH.out.reapr_per_base_bedgraph )

	TABIX_BGZIPTABIX( UCSC_BIGWIGTOBEDGRAPH.out.bedgraph )

	TABIX_BGZIPTABIX.out.gz_tbi.filter{ v -> v[1].toString()=~/.*_base.*gz$/ }.map{ meta, bedgraph, bedgraph_tbi -> [['id': meta['id_original']], bedgraph] }.set{bedgraph_filtered_ch}
	TABIX_BGZIPTABIX.out.gz_tbi.filter{ v -> v[2].toString()=~/.*_base.*gz\.tbi$/ }.map{ meta, bedgraph, bedgraph_tbi -> [['id': meta['id_original']], bedgraph_tbi] }.set{bedgraph_filtered_tbi_ch}

	bedgraph_filtered_ch.view{ "\n--------------\nTEMPLATE TO REPORT: $it \n\n" }
	TABIX_BGZIPTABIX.out.gz_tbi.filter{ v -> v[1].toString()=~/.*_depth.*gz$/ }.map{ meta, bedgraph, bedgraph_tbi -> [['id': meta['id_original']], bedgraph] }.set{bedgraph_filtered_depth_ch}
	TABIX_BGZIPTABIX.out.gz_tbi.filter{ v -> v[2].toString()=~/.*_depth.*gz\.tbi$/ }.map{ meta, bedgraph, bedgraph_tbi -> [['id': meta['id_original']], bedgraph_tbi] }.set{bedgraph_filtered_depth_tbi_ch}

	TABIX_BGZIPTABIX.out.gz_tbi.filter{ v -> v[1].toString()=~/.*_insert.*gz$/ }.map{ meta, bedgraph, bedgraph_tbi -> [['id': meta['id_original']], bedgraph] }.set{bedgraph_filtered_insert_ch}
	TABIX_BGZIPTABIX.out.gz_tbi.filter{ v -> v[2].toString()=~/.*_insert.*gz\.tbi$/ }.map{ meta, bedgraph, bedgraph_tbi -> [['id': meta['id_original']], bedgraph_tbi] }.set{bedgraph_filtered_insert_tbi_ch}

	TABIX_BGZIPTABIX.out.gz_tbi.filter{ v -> v[1].toString()=~/.*_kmer.*gz$/ }.map{ meta, bedgraph, bedgraph_tbi -> [['id': meta['id_original']], bedgraph] }.set{bedgraph_filtered_kmer_ch}
	TABIX_BGZIPTABIX.out.gz_tbi.filter{ v -> v[2].toString()=~/.*_kmer.*gz\.tbi$/ }.map{ meta, bedgraph, bedgraph_tbi -> [['id': meta['id_original']], bedgraph_tbi] }.set{bedgraph_filtered_kmer_tbi_ch}

	TABIX_BGZIPTABIX.out.gz_tbi.filter{ v -> v[1].toString()=~/.*_place.*gz$/ }.map{ meta, bedgraph, bedgraph_tbi -> [['id': meta['id_original']], bedgraph] }.set{bedgraph_filtered_place_ch}
	TABIX_BGZIPTABIX.out.gz_tbi.filter{ v -> v[2].toString()=~/.*_place.*gz\.tbi$/ }.map{ meta, bedgraph, bedgraph_tbi -> [['id': meta['id_original']], bedgraph_tbi] }.set{bedgraph_filtered_place_tbi_ch}

	REAPR.out.reapr_bam.map{ meta_complete, reapr_bam -> [['id': meta_complete['id']], reapr_bam]}.set{ reapr_bam_small_meta_ch }
	REAPR.out.reapr_bai.map{ meta_complete, reapr_bai -> [['id': meta_complete['id']], reapr_bai]}.set{ reapr_bai_small_meta_ch }

	TABIX_BGZIPTABIX_REAPR.out.gz_tbi.map{ meta, bedgraph, bedgraph_tbi -> [meta, bedgraph] }.set{bedgraph_filtered_reapr_per_base_ch}
	TABIX_BGZIPTABIX_REAPR.out.gz_tbi.map{ meta, bedgraph, bedgraph_tbi -> [meta, bedgraph_tbi] }.set{bedgraph_filtered_reapr_per_base_tbi_ch}

	bedgraph_filtered_reapr_per_base_ch.view{ "\n@@@@@@@@@@@@@@@@@@@@@@@@nREAPER BEDGRAPH: $it \n\n" }


	template_track_igv.join(reapr_dir_small_meta_ch).join(bedgraph_filtered_ch)
		.join(bedgraph_filtered_depth_ch).join(bedgraph_filtered_insert_ch)
		.join(bedgraph_filtered_kmer_ch).join(bedgraph_filtered_place_ch)
		.join(reapr_bam_small_meta_ch).join(reapr_bai_small_meta_ch)
		.join(bedgraph_filtered_reapr_per_base_ch)
		.multiMap{ meta_id, template, reapr_dir, ale_wig_base, ale_wig_depth, ale_wig_insert, ale_wig_kmer, ale_wig_place, reapr_bam, reapr_bai, reapr_score_per_base -> 
			meta_id: meta_id
			template: template
			reapr_dir: reapr_dir
			ale_wig_base: ale_wig_base
			ale_wig_depth: ale_wig_depth
			ale_wig_insert: ale_wig_insert
			ale_wig_kmer: ale_wig_kmer
			ale_wig_place: ale_wig_place
			reapr_bam: reapr_bam
			reapr_bai: reapr_bai
			reapr_score_per_base: reapr_score_per_base
		}.set{ template_to_report_ch }

	template_to_report_ch.meta_id.view{ "\n--------------\nTEMPLATE TO REPORT: $it \n\n" }
	

	PREPARE_REPORT_IGV ( template_to_report_ch.template, template_to_report_ch.meta_id, template_to_report_ch.reapr_bam, template_to_report_ch.reapr_bai,
	template_to_report_ch.ale_wig_base, template_to_report_ch.ale_wig_depth,
	template_to_report_ch.ale_wig_insert, template_to_report_ch.ale_wig_kmer,
	template_to_report_ch.ale_wig_place, template_to_report_ch.reapr_score_per_base )

	// PREPARE_REPORT_IGV.out.view{ "\n--------------\nPREPARE REPORT IGV: $it \n\n"}
	SUBSET_REAPR.out.reaper_bed_failure.map{ meta_complete, reapr_failure -> [['id': meta_complete['id']], reapr_failure]}.set{ reapr_failure_small_meta_ch }

	reapr_failure_small_meta_ch.join(PREPARE_REPORT_IGV.out.track_files)
	.join(input_chr_size_ch.asm)
	.join(CUSTOM_GETCHROMSIZES.out.fai)
	.join(PREPARE_REPORT_IGV.out.track_config)
	.join(bedgraph_filtered_tbi_ch)
	.join(bedgraph_filtered_depth_tbi_ch)
	.join(bedgraph_filtered_insert_tbi_ch)
	.join(bedgraph_filtered_kmer_tbi_ch)
	.join(bedgraph_filtered_place_tbi_ch)
	.join(bedgraph_filtered_reapr_per_base_tbi_ch)
	.multiMap{ meta_id, reapr_failure,
		 reapr_bam, reapr_bai, ale_wig_base, ale_wig_depth, ale_wig_insert, ale_wig_kmer, ale_wig_place, reaper_score_per_base, asm_fasta, asm_fai, track_config,
		 bedgraph_filtered_tbi_ch, bedgraph_filtered_depth_tbi_ch, bedgraph_filtered_insert_tbi_ch, bedgraph_filtered_kmer_tbi_ch,  bedgraph_filtered_place_tbi_ch, bedgraph_filtered_reapr_per_base_ch -> 
			in_tracks: [meta_id, reapr_failure, [ale_wig_base, ale_wig_depth, ale_wig_insert, ale_wig_kmer, ale_wig_place, reapr_bam, reapr_bai, reaper_score_per_base],[bedgraph_filtered_tbi_ch, bedgraph_filtered_depth_tbi_ch, bedgraph_filtered_insert_tbi_ch, bedgraph_filtered_kmer_tbi_ch,  bedgraph_filtered_place_tbi_ch, bedgraph_filtered_reapr_per_base_ch]]
			in_fasta: [meta_id, asm_fasta, asm_fai]
			in_config: [meta_id, track_config]
	}.set{ in_tracks_report_ch }


	IGVREPORTS ( in_tracks_report_ch.in_tracks, in_tracks_report_ch.in_fasta, in_tracks_report_ch.in_config )

	IGVREPORTS.out.report.view{ "\n--------------\n REPORT IGV: $it \n\n" }

	reapr_out_ch.mix( REAPR.out.reapr_summary ).set{ reapr_out_ch }

	igv_report_out_ch.mix( IGVREPORTS.out.report.map{ meta, report -> report } ).set{ igv_report_out_ch }

	// Execute QUAST without a reference and a GFF
	// QUAST ( reference_ch, [[], []], [[], []] );

    // if( params.aligner == 'bowtie' ){
	// 	inputs_reads_asm_ch.map{ meta_cor_reads, reads, meta_asm, asm -> [ asm ] }.set{ reference_bowtie_ch }
	// 	BOWTIE_BUILD ( reference_bowtie_ch )

	// 	inputs_reads_asm_ch.map{ meta_cor_reads, reads, meta_asm, asm -> [ meta_cor_reads, reads ] }.set{ reads_ch_cor_meta }
	// 	BOWTIE_ALIGN ( reads_ch_cor_meta, BOWTIE_BUILD.out.index )

	// 	aligned_reads_ch.mix( BOWTIE_ALIGN.out.bam )
    //         .set { aligned_reads_ch }
    // } else if ( params.aligner == 'bwa' ) {
	// 	inputs_reads_asm_ch.map{ meta_cor_reads, reads, meta_asm, asm -> [ meta_asm, asm ] }.set{ reference_bwa_ch }
    //     BWA_INDEX ( reference_bwa_ch )

	// 	inputs_reads_asm_ch.map{ meta_cor_reads, reads, meta_asm, asm -> [ meta_cor_reads, reads ] }.set{ reads_ch_cor_meta }
    //     BWA_MEM ( reads_ch_cor_meta, BWA_INDEX.out.index, true )

	// 	aligned_reads_ch.mix( BWA_MEM.out.bam )
    //         .set { aligned_reads_ch }
    // }
    // aligned_reads_ch.view()

    emit:
    ale = ale_out_ch   // queue channel: [ sample_id, file(bam_file) ]
	reapr = reapr_out_ch
	igv_report = igv_report_out_ch
	versions = REAPR.out.versions.mix(ALE.out.versions)
}