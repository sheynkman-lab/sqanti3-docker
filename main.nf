
if (!params.sample_gtf) exit 1, "Cannot find any file for parameter --sample_gtf: ${params.sample_gtf}"
ch_sample_gtf =  Channel.value(file(params.sample_gtf))

if (!params.best_orf) exit 1, "Cannot find any file for parameter --best_orf: ${params.best_orf}"
ch_best_orf =  Channel.value(file(params.best_orf))

if (!params.reference_gtf) exit 1, "Cannot find any file for parameter --reference_gtf: ${params.reference_gtf}"
ch_reference_gtf =  Channel.value(file(params.reference_gtf))

process rename_cds_to_exon{
    publishDir "${params.outdir}/${params.name}/rename/", mode: 'copy'
    tag "${params.name} ${reference_gtf} ${sample_gtf}"
    input:
        file(reference_gtf) from ch_reference_gtf
        file(sample_gtf) from ch_sample_gtf
    output:
        // file("*")
        file("${params.name}.cds_renamed_exon.gtf") into ch_sample_cds_renamed
        file("${params.name}.transcript_exons_only.gtf") into ch_sample_transcript_exon_only
        file("gencode.cds_renamed_exon.gtf") into ch_ref_cds_renamed
        file("gencode.transcript_exons_only.gtf") into ch_ref_transcript_exon_only

    script:
        """
        rename_cds_to_exon.py \
        --sample_gtf $sample_gtf \
        --sample_name ${params.name} \
        --reference_gtf $reference_gtf \
        --reference_name gencode
        """
}

process sqanti_protein{
    publishDir "${params.outdir}/${params.name}/sqanti_protein/", mode: 'copy'
    input:
        file(sample_exon) from ch_sample_transcript_exon_only
        file(sample_cds) from ch_sample_cds_renamed
        file(reference_exon) from ch_ref_transcript_exon_only
        file(reference_cds) from ch_ref_cds_renamed
        file(best_orf) from ch_best_orf
    output:
        file("${params.name}.sqanti_protein_classification.tsv") into ch_sqanti_protein_classification
    script:
    """
    sqanti3_protein.py \
    $sample_exon \
    $sample_cds \
    $best_orf \
    $reference_exon \
    $reference_cds \
    -d ./ \
    -p ${params.name}
    """
}


process protein_classification{
    publishDir "${params.outdir}/${params.name}/protein_classification/", mode: 'copy'
    input:
        file(sqanti_protein) from ch_sqanti_protein_classification
    output:
        file("${params.name}.protein_classification.tsv")
    script:
        """
        protein_classification.py \
        --sqanti_protein $sqanti_protein \
        --name ${params.name}
        """

}
