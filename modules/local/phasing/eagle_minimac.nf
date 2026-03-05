process EAGLE_MINIMAC {

    label 'phasing_imputation'
    tag "${chunkfile}"
    
    input:
    tuple val(chr), path(bcf), path(bcf_csi), val(start), val(end), val(phasing_status), path(chunkfile), path(m3vcf)
    path map_eagle
    path minimac_map
    val refpanel_build
    val minimac_window
    val minimac_min_ratio
    val min_r2
    val decay
    val diffThreshold
    val probThreshold
    val probThresholdS1
    val minRecombination

    output:
    tuple val(chr_cleaned), val(start), val(end), file("*.dose.vcf.gz"), file("*.info.gz"), file("*.empiricalDose.vcf.gz"), emit: em_imputed_chunks

    script:
    //define basename without ending (do not use simpleName due to X.*)
    def chunkfile_name = chunkfile.toString().replaceAll('.vcf.gz', '')
    // replace X.nonPAR etc with X for phasing
    def chr_cleaned = chr.startsWith('X.') ? 'X' : chr
    def chr_mapped = params.refpanel.build == 'hg38' ? 'chr' + chr_cleaned : chr_cleaned
    def phasing_start = start.toLong() - params.phasing.window
    phasing_start = phasing_start < 0 ? 1 : phasing_start
    def phasing_end = end.toLong() + params.phasing.window
    def used_threads = params.service.threads != -1 ? params.service.threads : task.cpus
    def chunkfile2=path(${chunkfile_name}.phased.vcf.gz)
    // minimac4
    map = minimac_map ? '--map ' + minimac_map : ''
    r2_filter = min_r2 != 0 ? '--min-r2 ' + min_r2 : ''
    diff_threshold = diffThreshold != -1 ? '--diff-threshold ' + diffThreshold : ''
    prob_threshold = probThreshold != -1 ? '--prob-threshold ' + probThreshold : ''
    prob_threshold_s1 = probThresholdS1 != -1 ? '--prob-threshold-s1 ' + probThresholdS1 : ''
    min_recom = minRecombination != -1 ? '--min-recom ' + minRecombination : ''

    """
    tabix $chunkfile

    eagle --vcfRef ${bcf} --vcfTarget ${chunkfile} --geneticMapFile ${map_eagle} --outPrefix ${chunkfile_name}.phased --chrom $chr_mapped --bpStart $phasing_start --bpEnd $phasing_end --allowRefAltSwap --vcfOutFormat z --keepMissingPloidyX --numThreads $used_threads

    tabix $chunkfile2

    minimac4 --region $chr_mapped:$start-$end --overlap $minimac_window --format GT,DS,GP,HDS --min-ratio $minimac_min_ratio --all-typed-sites --sites ${chunkfile_name}.info.gz --empirical-output ${chunkfile_name}.empiricalDose.vcf.gz --output ${chunkfile_name}.dose.vcf.gz --output-format vcf.gz --threads $used_threads --decay $decay --temp-prefix ./ $diff_threshold $prob_threshold $prob_threshold_s1 $min_recom $r2_filter $map ${m3vcf} ${chunkfile2}
    """
}
