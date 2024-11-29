/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_assemblyeval_pipeline'

// def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
// def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
// def summary_params = paramsSummaryMap(workflow)

// // Print parameter summary log to screen
// log.info logo + paramsSummaryLog(workflow) + citation

// WorkflowAssemblyeval.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
// include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { PREPARE_INPUT } from '../subworkflows/local/prepare_input'
include { READ_MAPPING } from '../subworkflows/local/read_mapping'
include { CONTIGUITY_ASM } from '../subworkflows/local/contiguity'
include { COMPLETENESS_ASM } from '../subworkflows/local/completeness'
include { CORRECTNESS_ASM } from '../subworkflows/local/correctness'
include { KMER_PROFILE } from '../subworkflows/local/kmerprofile'
include { CONTAMINATION_ASM } from '../subworkflows/local/contamination'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
// include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { PARSE_RESULTS } from '../modules/local/parse_results'
include { JINJA_MULTIQC } from '../modules/local/jinja_multiqc'
include { JINJA_CONTAMINANT } from '../modules/local/jinja_contaminant_multiqc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
// def multiqc_report = []

workflow ASSEMBLYEVAL {

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_asm = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    PREPARE_INPUT (
        Channel.fromPath(params.input)
    )
    
    ch_versions = ch_versions.mix(PREPARE_INPUT.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    illumina_ch = PREPARE_INPUT.out.illumina
    
    // illumina_ch.view{"READS CHANNEL: " + it}

    // PREPARE_INPUT.out.assemblies.view{"\n\nASSEMBLIES CHANNEL: " + it}

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        illumina_ch
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
     
    PREPARE_INPUT.out.assemblies.map{ meta, asm -> [meta, asm]}.view{ "OUT ASSEMBLIES: $it" }
    //
    // MODULE: Contamination
    //
    CONTAMINATION_ASM( illumina_ch, PREPARE_INPUT.out.assemblies )

    ch_asm = CONTAMINATION_ASM.out.asms.filter{ it != [] }

    ch_asm.view{ "OUT MIXED: $it" }

    Channel.fromPath("$projectDir/assets/template_contaminant_mqc.html").combine(CONTAMINATION_ASM.out.fcs).set{ to_jinja_contaminant_ch }

    JINJA_CONTAMINANT (
        to_jinja_contaminant_ch.map{ tp, meta, fcs_report -> tp },
        to_jinja_contaminant_ch.map{ tp, meta, fcs_report -> meta.get("id") },
        to_jinja_contaminant_ch.map{ tp, meta, fcs_report -> fcs_report }
    )


    JINJA_CONTAMINANT.out.rendered.view{ "================== OUT JINJA CONTAMINANT: $it"}

    ch_multiqc_files = ch_multiqc_files.mix(JINJA_CONTAMINANT.out.rendered)

    //
    // MODULE: Completeness Metrics
    //
    COMPLETENESS_ASM (
        illumina_ch, ch_asm
    )

    // 
    // MODULE: Contiguity Metrics
    // 
    CONTIGUITY_ASM (
        ch_asm
    )

    // //
    // // MODULE: Correctness Metrics
    // //
    // // SUB-MODULE: Mapping reads
    // //
    READ_MAPPING (
        illumina_ch, ch_asm
    )

    // // READ_MAPPING.out.bam.view{ "*** ETAPA 2 - BAM: $it"}
    // // READ_MAPPING.out.bai.view{ "*** ETAPA 2 - BAI: $it"}

    
    CORRECTNESS_ASM (
       READ_MAPPING.out.joined_asm_bam, READ_MAPPING.out.asm, READ_MAPPING.out.bam, READ_MAPPING.out.bai
    )

    KMER_PROFILE (
        illumina_ch, ch_asm
    )

    // PREPARE_INPUT.out.assemblies.view{ "ASSEMBLIES TO KMER PROFILE: $it" }
    // // parse_results_to_table(args.genomes_ids, args.ale_res, args.reapr_res, args.busco_re_summary, args.quast_res, args.file_out)
    
    out_asm_ch = CORRECTNESS_ASM.out.reapr.map{ meta, reapr -> meta['id'] }.collect( sort: true ).view{ "ASSEMBLIES: $it" }
    out_ale_ch = CORRECTNESS_ASM.out.ale.toSortedList( { a, b -> a[0]['id'] <=> b[0]['id'] } ).flatMap().collect{ it[1] }
    out_reapr_ch = CORRECTNESS_ASM.out.reapr.toSortedList( { a, b -> a[0]['id'] <=> b[0]['id'] } ).flatMap().collect{ it[1] }
    out_busco_re_summary_ch = COMPLETENESS_ASM.out.busco_short_summaries_txt.toSortedList( { a, b -> a[0]['id'] <=> b[0]['id'] } ).flatMap().collect{ it[1] }
    out_quast_ch = CONTIGUITY_ASM.out.quast_tsv.toSortedList( { a, b -> a[0]['id'] <=> b[0]['id'] } ).flatMap().collect{ it[1] }
    out_merfin_qv_ch = KMER_PROFILE.out.merfin_logs.toSortedList( { a, b -> a[0]['id'] <=> b[0]['id'] } ).flatMap().collect{ it[1] }
    out_merfin_completeness_ch = KMER_PROFILE.out.merfin_completeness.toSortedList( { a, b -> a[0]['id'] <=> b[0]['id'] } ).flatMap().collect{ it[1] }

    // // out_table_ch = Channel.fromPath("out_table.csv")
    
    // COMPLETENESS_ASM.out.busco_short_summaries_txt.view{ "REAPR: $it" }
    
    PARSE_RESULTS ( out_asm_ch, out_ale_ch, out_reapr_ch, out_busco_re_summary_ch, out_quast_ch, out_merfin_qv_ch, out_merfin_completeness_ch  )
    // PARSE_RESULTS.out.res.view{ "\n\nRESULTS ARE IN: $it"}

    // // // CORRECTNESS_ASM.out.reapr.collect( {it[1]}, sort: {it.getName()} ).set{ new_collect }
    

    // // ch_versions = ch_versions.mix(FASTQC.out.versions.first()).mix(CONTIGUITY_ASM.out.versions.first())


    // CUSTOM_DUMPSOFTWAREVERSIONS (
    //     ch_versions.unique().collectFile(name: 'collated_versions.yml')
    // )

    ch_versions = ch_versions.mix(KMER_PROFILE.out.versions).mix(CORRECTNESS_ASM.out.versions).mix(READ_MAPPING.out.versions).mix(CONTIGUITY_ASM.out.versions).mix(COMPLETENESS_ASM.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    ch_multiqc_custom_config = ch_multiqc_custom_config.mix(Channel.fromPath("$projectDir/assets/multiqc_kmer_config.yaml", checkIfExists: true))
    ch_multiqc_custom_config = ch_multiqc_custom_config.mix(Channel.fromPath("$projectDir/assets/multiqc_contaminant_config.yaml", checkIfExists: true))

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )
    ch_multiqc_files = ch_multiqc_files.mix(PARSE_RESULTS.out.res.map{ asm_names, table -> table})
    // ch_multiqc_files = ch_multiqc_files.mix(Channel.fromPath("$projectDir/assets/"))


    KMER_PROFILE.out.merqury_cn.map{ meta, cn_spectra -> [['to_join': meta.get("to_join")], meta, cn_spectra] }.set{ to_jinja_merqury_ch }
    view{ "MERQURY CN: $it" }

    KMER_PROFILE.out.genomescope2_res.map{ meta, gs_ln, gs_fitted -> [['to_join': meta.get('id')], meta, gs_ln, gs_fitted] }.set{ to_jinja_genscope_ch }

    to_jinja_genscope_ch.combine(to_jinja_merqury_ch, by: [0]).map{ meta_join, meta, gs_ln , gs_fitted, meta_mq, cn_spectra -> [meta_mq, gs_ln, cn_spectra] }.set{ to_jinja_samp_ch }

    Channel.fromPath("$projectDir/assets/template_mqc.html").combine(to_jinja_samp_ch).set{ to_jinja_ch }

    // to_jinja_ch.view{ "TO JINJA: $it" }

    JINJA_MULTIQC ( to_jinja_ch.map{ tp, meta, gs_ln, cn_spectra -> tp },
    to_jinja_ch.map{ tp, meta, gs_ln, cn_spectra -> meta.get("id") },
    to_jinja_ch.map{ tp, meta, gs_ln, cn_spectra -> gs_ln },
    to_jinja_ch.map{ tp, meta, gs_ln, cn_spectra  -> cn_spectra } )


    
     
    // ch_multiqc_files = ch_multiqc_files.mix(Channel.fromPath("$projectDir/assets/html_with_comment_mqc.html"))
    ch_multiqc_files = ch_multiqc_files.mix(JINJA_MULTIQC.out.rendered)
    
    ch_multiqc_logo = ch_multiqc_logo.mix(Channel.fromPath("$projectDir/assets/nf-core-assemblyeval_logo_light.png"))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    // //
    // // MODULE: MultiQC
    // //
    // workflow_summary    = WorkflowAssemblyeval.paramsSummaryMultiqc(workflow, summary_params)
    // ch_workflow_summary = Channel.value(workflow_summary)

    // methods_description    = WorkflowAssemblyeval.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    // ch_methods_description = Channel.value(methods_description)

    // ch_multiqc_files = Channel.empty()
    // ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    // ch_multiqc_files = ch_multiqc_files.mix(CONTIGUITY_ASM.out.quast.collect{it[1]}.ifEmpty([]))
    // if ( params.busco ) {
    //     ch_multiqc_files = ch_multiqc_files.mix(COMPLETENESS_ASM.out.busco_short_summaries_txt.collect{it[1]}.ifEmpty([]))
    // }
    // ch_multiqc_files = ch_multiqc_files.mix( PARSE_RESULTS.out.res.collect{it[1]}.ifEmpty([]))


    // // ch_multiqc_files.collect().view{ "QUAST_TSV: $it" }

    // ch_multiqc_custom_config.mix( Channel.fromPath("./assets/section_name_with_slash.yml")).set{ ch_multiqc_custom_config }

    // MULTIQC (
    //     ch_multiqc_files.collect(),
    //     ch_multiqc_config.toList(),
    //     ch_multiqc_custom_config.toList(),
    //     ch_multiqc_logo.toList()
    // )
    multiqc_report = MULTIQC.out.report.toList()
    emit:
    multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {

    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()
        
    // sendMail(to: 'theodoro.biotec@gmail.com', from: 'theodoro.biotec@gmail.com', subject: 'My pipeline execution', body: msg)
}

// workflow.onComplete {
//     if (params.email || params.email_on_fail) {
//         NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
//     }
//     NfcoreTemplate.dump_parameters(workflow, params)
//     NfcoreTemplate.summary(workflow, params, log)
//     if (params.hook_url) {
//         NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
//     }
// }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
