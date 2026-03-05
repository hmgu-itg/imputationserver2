include { EAGLE_MINIMAC } from '../modules/local/phasing/eagle_minimac'

workflow PHASING_IMPUTATION {
    take: 
    metafiles_ch 
    main:

    chromosomes = Channel.of(1..22, 'X.nonPAR', 'X.PAR1', 'X.PAR2', 'MT')
    if (params.phasing.engine == 'eagle' && params.refpanel.refEagle != null) {
        phasing_map_ch = file(params.refpanel.mapEagle, checkIfExists: true)
        minimac_map = []
        if (params.refpanel.mapMinimac != null) {
            minimac_map = file(params.refpanel.mapMinimac, checkIfExists: true)
        }
        phasing_reference_ch = chromosomes
            .map {
                it -> 
                    def eagle_file = file(PatternUtil.parse(params.refpanel.refEagle, [chr: it]))
                    def eagle_file_index = file(PatternUtil.parse(params.refpanel.refEagle + ".csi", [chr: it]))
                    if(!eagle_file.exists() || !eagle_file_index.exists()){
                        return null;
                    }
                    return tuple(it.toString(),eagle_file,eagle_file_index)
            }
	minimac_m3vcf_ch = chromosomes
        .map {
            it ->
                def genotypes_file = file(PatternUtil.parse(params.refpanel.genotypes, [chr: it]))
                    if(!genotypes_file.exists()){
                        return null;
                    }
                return tuple(it.toString(),genotypes_file);
        }
        eagle_bcf_metafiles_ch = phasing_reference_ch.combine(metafiles_ch, by: 0)
        em_metafiles_ch = eagle_bcf_metafiles_ch.combine(minimac_m3vcf_ch, by: 0)
        EAGLE_MINIMAC ( em_metafiles_ch, phasing_map_ch,minimac_map,
        params.refpanel.build,
        params.imputation.window,
        params.imputation.minimac_min_ratio,
        params.imputation.min_r2,
        params.imputation.decay,
        params.imputation.diff_threshold,
        params.imputation.prob_threshold,
        params.imputation.prob_threshold_s1,
        params.imputation.min_recom
	)

    emit:
    imputed_chunks = EAGLE_MINIMAC.out.em_imputed_chunks

}
