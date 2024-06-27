include { CAT_FASTQ } from '../../modules/nf-core/cat/fastq/main'
include { MERYL_COUNT as MERYL_COUNT_READS_01;
          MERYL_COUNT as MERYL_COUNT_GENOME_01  } from '../../modules/nf-core/meryl/count/main'
include { MERYL_GREATER_THAN } from '../../modules/local/meryl_greater'
include { MERYL_HISTOGRAM as MERYL_HISTOGRAM_READS_PRE;
          MERYL_HISTOGRAM as MERYL_HISTOGRAM_GENOME_PRE } from '../../modules/nf-core/meryl/histogram/main'
include { GENOMESCOPE2 as GENOMESCOPE2_PRE } from '../../modules/nf-core/genomescope2/main'
include { MERFIN_COMPLETENESS } from '../../modules/local/merfin/completeness/main'
include { MERFIN_HIST as MERFIN_HIST_EVALUATE_POLISH;
          MERFIN_HIST as MERFIN_HIST_EVALUATE_PRIMARY_ASSEMBLY } from '../../modules/nf-core/merfin/hist/main'

process GET_PEAK {
    input:
    tuple val(meta), path(genomescope2model)

    output:
    tuple val(meta), path("*_peak.txt"), emit: peak
// cat ${genomescope2model} | grep "^kmercov" | python -c 'import sys; print(float([x.split() for x in sys.stdin.readlines()][0][1]))' > file

    // peak="\$(cat ${genomescope2model} | grep "^kmercov" | python -c 'import sys; print(float([x.split() for x in sys.stdin.readlines()][0][1]))')"
    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat ${genomescope2model} | grep "^kmercov" | python -c 'import sys; print(float([x.split() for x in sys.stdin.readlines()][0][1]))' > ${prefix}_peak.txt
    """
}

workflow KMER_PROFILE {

	take:
    reads        // queue channel; [ sample_id, [ file(read1), file(read2) ] ]
    reference    // file( "path/to/reference" )

    main:

	// (0) Prepare nf empty variables with will receive results
    Channel.empty()
        .set { aligned_reads_ch }

	Channel.empty()
        .set { res_completeness_ch }

	Channel.empty()
        .set { res_logs_ch }
	
    Channel.empty()
        .set { res_versions_ch }
    
	// (1) Generate k-mers count distribution

	// (1.1) Count operation: count the occurrences of canonical k-mers in high-accuracy reads
	
	// Create a channel where paired-end data is mentioned as single-end
    // This is necessary in order to concatenate paired-end in a single file by cat/fastq nf-core module
	
    reads_ch_single = reads.map { meta, reads -> {
		def meta_new = meta.clone();
		meta_new.put('single_end', true);
		[meta_new, reads] }
	}
	// reads_ch_single.view{ "++++++++++++ READS: $it" }
            
    CAT_FASTQ (
        reads_ch_single
    )
	
    CAT_FASTQ.out.reads.multiMap{ meta, reads -> 
        outs: [meta, reads]
        kmer: meta['kmer_size']
    }.set{ meryl_ch }

	MERYL_COUNT_READS_01 (
        meryl_ch.outs, meryl_ch.kmer
    )

	// (1.2) Rename output as meryldb - SKIPPED BECAUSE WE JOIN THE READS?
	// (1.3) Operation Union-sum: return k-mers that occur in any input, set the count to the sum of the counts - SKIPPED BECAUSE WE JOIN THE READS?
	// (1.4) Rename it as merged_meryldb - SKIPPED BECAUSE WE JOIN THE READS?

	// (1.5-ROD) Keep only k-mers with frequency > 1
	MERYL_GREATER_THAN (
        MERYL_COUNT_READS_01.out.meryl_db
    )

	// (1.6) Generate histogram 
	MERYL_HISTOGRAM_READS_PRE (
        // MERYL_COUNT_READS_01.out.meryl_db.mix(MERYL_COUNT_GENOME_01.out.meryl_db)
        // MERYL_GREATER_THAN.out.meryl_db
		MERYL_COUNT_READS_01.out.meryl_db, meryl_ch.kmer
    )

	// (1.7) Genome profiling with GenomeScope2
	GENOMESCOPE2_PRE (
        MERYL_HISTOGRAM_READS_PRE.out.hist
    )

    peak_out = GET_PEAK (
        GENOMESCOPE2_PRE.out.model
    )

    peak_ch_val = peak_out.map{ meta, peak -> peak.text.trim() }.first()
    // peak_ch_val_test = peak_out.map{ it.text.trim() }.first()

	// (1.9) Genome assembly k-mer count
	reference.map{ meta, assemblies -> {
		def meta_new = meta.clone();
		meta_new.put("id", meta.get("id")+"_"+meta.get("build"));
		[ ['id_to_join': meta.get("id")], meta_new, assemblies['pri_asm'] ]
	 } }.set{ reference_ch }

    // reference_ch.view{"ASSEMBLIES MEFIN: $it"}

	MERYL_COUNT_GENOME_01 (
        reference_ch.map{ meta_join, meta, assembly -> [ meta, assembly] }, meryl_ch.kmer
    )

    // MERYL_GREATER_THAN.out.meryl_db.combine(GENOMESCOPE2_PRE.out.lookup_table).combine(peak_ch_val).set{ to_merfin_ch }
    MERYL_GREATER_THAN.out.meryl_db.combine(GENOMESCOPE2_PRE.out.lookup_table, by: [0]).combine(peak_out, by: [0]).map{ meta, readmer, lt, peak ->
         [['id_to_join': meta.get("id")], meta, readmer, meta, lt, peak.text.trim()] }.set{ to_merfin_ch }

    // to_merfin_ch_test.view{ "TO MERFIN: $it" }
    
    // to_merfin_ch.map{ meta_reads_db, meryl_db, meta_lt, lt, peak -> 
    // [ ['id_to_join': meta_reads_db.get("id")], meta_reads_db, meryl_db, meta_lt, lt, peak] }.set{ to_merfin_ch }

    // reference_ch.combine(to_merfin_ch_test, by: [0]).map{ meta_join, meta, asm, meta_reads_db, meryl_db, meta_lt, lt, peak -> 
    //     [[meta, asm],[meta_reads_db, meryl_db],lt,[],peak]  }.set{ to_test_merfin_ch }

    // to_test_merfin_ch.view{ "TO MERFIN: $it" }

    reference_ch.combine(to_merfin_ch, by: [0]).multiMap{ meta_join, meta, asm, meta_reads_db, meryl_db, meta_lt, lt, peak -> 
        assembly: [meta, asm]
        readmers: [meta_reads_db, meryl_db]
        lookup_table: lt
        seqmers: []
        peak: peak  }.set{ to_merfin_ch }


    // to_test_merfin_ch.view{ "&&&&&&&&&&&&&&&& TO MERFIN ASM: $it" }
    // to_merfin_ch.map{ }.readmers.view{ "&&&&&&&&&&&&&&&& TO MERFIN READMER: $it\n\n" }


	// (1.8) K-mer based evaluation with Merfin (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9745813/)
    // MERYL_COUNT_GENOME_01.out.meryl_db,
	MERFIN_COMPLETENESS (
        to_merfin_ch.assembly,
        to_merfin_ch.readmers,
        to_merfin_ch.lookup_table,
        to_merfin_ch.seqmers,
        to_merfin_ch.peak
    )

    // MERFIN_COMPLETENESS.out.completeness.view{ "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& MERFIN - COMPLETENESS: $it" }

	// res_completeness_ch.mix(MERFIN_COMPLETENESS.out.completeness).set{ res_completeness_ch }	
    res_completeness_ch = MERFIN_COMPLETENESS.out.completeness
    // res_completeness_ch.view{ "MERFIN - COMPLETENESS: $it" }

    res_versions_ch.mix(MERFIN_COMPLETENESS.out.versions).set{ res_versions_ch }

    MERFIN_HIST_EVALUATE_PRIMARY_ASSEMBLY (
        to_merfin_ch.assembly,
        to_merfin_ch.readmers,
        to_merfin_ch.lookup_table,
        to_merfin_ch.seqmers,
        to_merfin_ch.peak
    )

    // MERFIN_HIST_EVALUATE_PRIMARY_ASSEMBLY.out.log_stderr.view{ "============================ MERFIN - MERYL DB: $it" }


    MERFIN_HIST_EVALUATE_PRIMARY_ASSEMBLY.out.hist.merge(MERFIN_HIST_EVALUATE_PRIMARY_ASSEMBLY.out.log_stderr).set{ out_merfin_hist_log_ch }


	




    

	// reference.map{ meta, assemblies -> {
	// 	def meta_new = meta.clone();
	// 	meta_new.put("id", meta.get("id")+"_"+meta.get("build"));
	// 	[ meta_new, assemblies['pri_asm'] ]
	//  } }.set{ reference_ch }
	

	// reference_ch.combine(reads).map{ meta_asm, asm, meta_reads, reads -> {
	// 	def meta_cor_reads = meta_reads.clone();
	// 	meta_cor_reads.put("id", meta_asm.get("id"));
	// 	[ meta_cor_reads, reads, meta_asm, asm ]
	// }
	// }.set{ inputs_reads_asm_ch }



    emit:
	res = aligned_reads_ch
	merfin_completeness = res_completeness_ch
	merfin_logs = out_merfin_hist_log_ch.map{ meta, hist, log -> [meta, log] }
    // quast = QUAST.out.results   // queue channel: [ sample_id, file(bam_file) ]
	// quast_tsv = QUAST.out.tsv
	// busco = BUSCO.out.busco_dir
	// busco_short_summaries_txt = BUSCO.out.short_summaries_txt
	versions = res_versions_ch
}