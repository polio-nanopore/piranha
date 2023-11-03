process run_piranha {

    publishDir "${params.out_dir}", mode: 'copy'

    container "${params.wf.container}:${workflow.manifest.version}"

    input:
        path barcodes_csv
        path run_dir

    output:
        path "piranha_report.html"
        path "barcode_reports"
        path "detailed_run_report.csv"
        path "published_data"

    script:
    extra = ""
    if ( params.config )
        extra += " --config ${params.config}"
    if ( params.output_intermediates )
        extra += " --no-temp"
    if ( params.min_map_quality )
        extra += " --min-map-quality ${params.min_map_quality}"
    if ( params.min_read_length )
        extra += " --min-read-length ${params.min_read_length}"
    if ( params.max_read_length )
        extra += " --max-read-length ${params.max_read_length}"
    if ( params.min_read_depth )
        extra += " --min-read-depth ${params.min_read_depth}"
    if ( params.min_read_pcent )
        extra += " --min-read-pcent ${params.min_read_pcent}"
    if ( params.primer_length )
        extra += " --primer-length ${params.primer_length}"
    if ( params.run_phylo )
        extra += " --run-phylo"
    if ( params.supplementary_sequences )
        extra += " --supplementary-sequences ${params.supplementary_sequences}"
     if ( params.supplementary_metadata )
            extra += " --supplementary-metadata ${params.supplementary_metadata}"
    """
    piranha -b ${barcodes_csv} -i ${run_dir} -o piranha_output --tempdir piranha_tmp -t ${task.cpus} ${extra}
    mv piranha_output/* .
    rm -rf piranha_output
    mv report.html piranha_report.html
    """

}

workflow {
    barcodes_csv = Channel.of(file("${params.barcodes_csv}", type: "file", checkIfExists:true))
    run_dir = Channel.of(file("${params.run_dir}", type: "dir", checkIfExists:true))

    run_piranha(barcodes_csv, run_dir)
}
