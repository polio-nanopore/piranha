process run_piranha {

    publishDir "${params.out_dir}", mode: 'copy'

    container "${params.wf.container}"

    input:
        path barcodes_csv
        path run_dir

    output:
        path "piranha_output"

    script:
    """
    piranha -b ${barcodes_csv} -i ${run_dir} -o piranha_output --tempdir piranha_tmp -t ${task.cpus}
    """

}

workflow {
    barcodes_csv = Channel.of(file("${params.barcodes_csv}", type: "file", checkIfExists:true))
    run_dir = Channel.of(file("${params.run_dir}", type: "dir", checkIfExists:true))

    run_piranha(barcodes_csv, run_dir)
}