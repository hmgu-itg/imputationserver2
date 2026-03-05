include { EAGLE } from '../modules/local/phasing/eagle'
include { BEAGLE } from '../modules/local/phasing/beagle'
include { MINIMAC4 } from '../modules/local/imputation/minimac4'

workflow PHASING_IMPUTATION {
    take: 
    metafiles_ch 
    main:

    chromosomes = Channel.of(1..22, 'X.nonPAR', 'X.PAR1', 'X.PAR2', 'MT')
    if (params.phasing.engine == 'eagle' && params.refpanel.refEagle != null) {
        phasing_map_ch = file(params.refpanel.mapEagle, checkIfExists: true)
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
        eagle_bcf_metafiles_ch = phasing_reference_ch.combine(metafiles_ch, by: 0)
    // imputation
    if (params.refpanel.mapMinimac == null) {
        minimac_map = []
    } else {
        minimac_map = file(params.refpanel.mapMinimac, checkIfExists: true)
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
    phased_m3vcf_ch = phased_ch.combine(minimac_m3vcf_ch, by: 0)
//        EAGLE ( eagle_bcf_metafiles_ch, phasing_map_ch )
//        phased_ch = EAGLE.out.eagle_phased_ch


    MINIMAC4 (
//        phased_m3vcf_ch,
EAGLE ( eagle_bcf_metafiles_ch, phasing_map_ch ).out.eagle_phased_ch.combine(minimac_m3vcf_ch, by: 0),
        minimac_map,
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

    imputed_chunks = MINIMAC4.out.imputed_chunks
    }

    // unmodified: does not include imputation
    // if (params.phasing.engine == 'beagle' && params.refpanel.refBeagle != null) {
    //     phasing_reference_ch = chromosomes
    //         .map {
    //             it -> 
    //                 def beagle_file = file(PatternUtil.parse(params.refpanel.refBeagle, [chr: it]))
    //                 if(!beagle_file.exists()){
    //                     return null;
    //                 }
    //                 return tuple(it.toString(),beagle_file)
    //         }
    //     phasing_map_ch = chromosomes
    //         .map {
    //             it ->
    //                 def beagle_map_file = file(PatternUtil.parse(params.refpanel.mapBeagle, [chr: it]))
    //                 if(!beagle_map_file.exists()){
    //                     return null;
    //                 }
    //                 return tuple(it.toString(),beagle_map_file)
    //         }
    //     beagle_bcf_metafiles_ch = phasing_reference_ch.combine(metafiles_ch, by: 0)
    //     //combine with map since also split by chromsome
    //     beagle_bcf_metafiles_map_ch = beagle_bcf_metafiles_ch.combine(phasing_map_ch, by: 0)
    //     BEAGLE ( beagle_bcf_metafiles_map_ch )
    //     phased_ch = BEAGLE.out.beagle_phased_ch
    // }

    emit:
    imputed_chunks

}
