#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    xtelin / Delly v1.7.2 BND detection Nextflow pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Sample  : 170TFW (SLC25A13 IVS16ins3kb breakpoint analysis)
    Tool    : Delly v1.7.2 (cross-chr BND + GT)
    Ref     : GRCh37/hg19 + hs37d5 + NC_007605
    Env     : docker/arm64 (xtelin/delly-bnd:v1.7.2)
----------------------------------------------------------------------------------------
    Usage:
      nextflow run main.nf -profile docker \
        --bam 170TFW.deduped.bam \
        --ref /mnt/storage/bioinformatics-data/kian/hg19/hs37d5.fa \
        --outdir ./nf_results

      # Batch (glob)
      nextflow run main.nf -profile docker \
        --bam '/path/to/samples/*.deduped.bam' \
        --ref /path/to/hs37d5.fa \
        --outdir ./nf_results
----------------------------------------------------------------------------------------
    Pipeline steps:
      1. ENSURE_BAM_INDEX   - Builds BAM .bai if missing
      2. ENSURE_FASTA_INDEX - Builds FASTA .fai, validates hs37d5/NC_007605 contigs
      3. DELLY_CALL         - delly call -t BND (whole-genome cross-chr BND + GT)
      4. EXTRACT_TARGET     - bcftools filter for chr6<->chr7 BNDs + target summary
      5. SUMMARY            - Per-sample final report
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INPUT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.bam    = null
params.ref    = null
params.outdir = 'nf_results'
params.extract_script = 'scripts/extract_bnd.py'

if (!params.bam) { error "Must specify --bam <file.bam> or --bam '*.bam'" }
if (!params.ref) { error "Must specify --ref <GRCh37_with_decoys.fa>" }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def bam_files = Channel.fromPath(params.bam, checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PROCESS: ENSURE_BAM_INDEX
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
process ENSURE_BAM_INDEX {
    tag { sample_name }
    container 'xtelin/delly-bnd:v1.7.2'

    input:
    tuple val(sample_name), path(bam)

    output:
    tuple val(sample_name), path(bam), path(bai)

    script:
    bai = bam.baseName + '.bai'
    """
    samtools index ${bam}
    if [ ! -f ${bai} ]; then
        echo "ERROR: BAM index failed: ${bai}"
        exit 1
    fi
    echo "BAM index OK: ${bai}"
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PROCESS: ENSURE_FASTA_INDEX
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
process ENSURE_FASTA_INDEX {
    container 'xtelin/delly-bnd:v1.7.2'

    input:
    path(fasta)

    output:
    tuple path(fasta), path(fai)

    script:
    fai = fasta.name + '.fai'
    """
    samtools faidx ${fasta}

    grep -q "^hs37d5"    ${fai} || echo "WARNING: missing hs37d5 contig"
    grep -q "^NC_007605" ${fai} || echo "WARNING: missing NC_007605 contig"
    echo "FASTA contigs: \\$(wc -l < ${fai})"
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PROCESS: DELLY_CALL
    Delly v1.7.2 cross-chr BND detection + genotyping
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
process DELLY_CALL {
    tag { sample_name }
    container 'xtelin/delly-bnd:v1.7.2'
    cpus   4
    memory '8G'

    input:
    tuple val(sample_name), path(bam), path(bai)
    tuple path(fasta), path(fai)

    output:
    tuple val(sample_name), path("${sample_name}.bcf"), path("${sample_name}.vcf")

    script:
    """
    echo "============================================"
    echo " Delly v1.7.2 cross-chr BND call"
    echo " Sample : ${sample_name}"
    echo " BAM    : ${bam}"
    echo " FASTA  : ${fasta}"
    echo "============================================"

    delly call -t BND -g ${fasta} -o ${sample_name}.bcf ${bam}

    bcftools view -O v -o ${sample_name}.vcf ${sample_name}.bcf
    bcftools index ${sample_name}.bcf

    echo ""
    echo "[Delly stats]"
    total=\\$(bcftools view ${sample_name}.vcf 2>/dev/null | grep -v '^#' | wc -l)
    pass=\\$(bcftools view -f 'PASS' ${sample_name}.vcf 2>/dev/null | grep -v '^#' | wc -l)
    echo "  total BND    : \\${total}"
    echo "  FILTER=PASS  : \\${pass}"
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PROCESS: EXTRACT_TARGET
    chr6<->chr7 BND extraction (SLC25A13 IVS16ins3kb breakpoint)
    Calls extract_bnd.py  outputs .bnd.json (6 fields, gene-aware)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
process EXTRACT_TARGET {
    tag { sample_name }
    container 'xtelin/delly-bnd:v1.7.2'

    input:
    tuple val(sample_name), path(bcf), path(vcf)
    path extract_script

    output:
    tuple val(sample_name), path("${sample_name}.target_bnd.vcf"),
                           path("${sample_name}.bnd.json")

    script:
    """
    # Delly  VCF  python   JSON 
    python3 ${extract_script} -t 6,7 --pass-only ${bcf}

    #  VCF chr6<->chr7 PASS BND
    bcftools view -T 6,7 -f 'PASS' ${bcf} 2>/dev/null \\
        | awk '\\$3=="#CHROM" || ((\\$1=="6" && /CHR2=7/) || (\\$1=="7" && /CHR2=6/))' \\
        > ${sample_name}.target_bnd.vcf

    echo ""
    echo "[JSON output]"
    cat ${sample_name}.bnd.json
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PROCESS: SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
process SUMMARY {
    container 'xtelin/delly-bnd:v1.7.2'

    input:
    tuple val(sample_name), path(vcf), path(json)

    output:
    path("pipeline_summary.txt")

    script:
    """
    {
        echo "=============================================="
        echo "  Delly BND Pipeline Summary"
        echo "=============================================="
        echo "Date    : \\$(date)"
        echo "Sample  : ${sample_name}"
        echo "----------------------------------------------"
        echo "  BND JSON (gene-aware, 6 fields):"
        cat ${json}
        echo "----------------------------------------------"
    } > pipeline_summary.txt
    cat pipeline_summary.txt
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow {
    // 1. Input channels (glob supported)
    bam_channel = bam_files.map { file ->
        def name = file.baseName.replaceAll(/\.deduped$/, '').replaceAll(/\.bam$/, '')
        [name, file]
    }

    // 2. FASTA index (single, shared)
    fasta_ch = ENSURE_FASTA_INDEX( file(params.ref) )

    // 3. BAM index (parallel per sample)
    indexed_bams = ENSURE_BAM_INDEX(bam_channel)

    // 4. Delly call (parallel per sample)
    delly_out = DELLY_CALL( indexed_bams, fasta_ch )

    // 5. Extract target BND    calls extract_bnd.py  .bnd.json
    target_out = EXTRACT_TARGET( delly_out, file(params.extract_script) )

    // 6. Final per-sample summary
    SUMMARY( target_out )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       PUBLISH RESULTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    delly_out.publishDir  "${params.outdir}/delly",  mode: 'copy', pattern: '*.{bcf,vcf,csi}'
    // target_out = tuple (sample_name, target_bnd.vcf, bnd.json)
    target_out.publishDir "${params.outdir}/target", mode: 'copy', pattern: '*.{vcf,json}'
    fasta_ch.publishDir   "${params.outdir}/ref",    mode: 'copy'
}
