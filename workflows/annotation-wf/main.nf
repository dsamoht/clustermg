include { BEDTOOLS                          } from '../../modules/bedtools' 
include { CHECKM                            } from '../../modules/checkm'
include { COLLECT                           } from '../../modules/collect'
include { DASTOOL                           } from '../../modules/dastool'
include { DASTOOL_CONTIG2BIN as METABAT_C2B } from '../../modules/dastool_contig2bin'
include { DASTOOL_CONTIG2BIN as MAXBIN_C2B  } from '../../modules/dastool_contig2bin'
include { FEATURECOUNTS                     } from '../../modules/featurecounts'
include { FEATURECOUNTS_SUMMARY             } from '../../modules/featurecounts_summary'
include { GTDBTK                            } from '../../modules/gtdbtk'
include { MAXBIN                            } from '../../modules/maxbin'
include { MAXBIN_ADJUST_EXT                 } from '../../modules/maxbin_adjust_ext'
include { METABAT                           } from '../../modules/metabat'
include { PRODIGAL                          } from '../../modules/prodigal'
include { SEQKIT                            } from '../../modules/seqkit'
include { TIARA                             } from '../../modules/tiara'
include { HMMER                             } from '../../modules/hmmer'
include { HMMER_SUMMARY                     } from '../../modules/hmmer_summary'


workflow ANNOTATION_WF {

    take:
    assembly
    sorted_bam
    
    main:
    PRODIGAL(assembly)

    FEATURECOUNTS(PRODIGAL.out.genesGff, sorted_bam)
    FEATURECOUNTS_SUMMARY(FEATURECOUNTS.out.counts)

    METABAT(assembly, sorted_bam)
    MAXBIN(assembly, METABAT.out.metabatDepth)
    MAXBIN_ADJUST_EXT(MAXBIN.out.maxbinBins)
    METABAT_C2B(METABAT.out.metabatBins, "metabat")
    MAXBIN_C2B(MAXBIN_ADJUST_EXT.out.renamed_maxbinBins, "maxbin")
        
    contig2bin_ch = METABAT_C2B.out.contigs2bins.
        mix(MAXBIN_C2B.out.contigs2bins).
        collect()

    DASTOOL(assembly, contig2bin_ch)
    TIARA(assembly)
    SEQKIT(DASTOOL.out.dasBins)

    if(params.profilePfam != '' || params.profileKegg != '') {
        profiles = Channel.of(["pfam", params.profilePfam], ["kegg", params.profileKegg])
        HMMER(genes = params.genesFaa, profiles.filter{ it.count('') == 0 })
        HMMER_SUMMARY(hmmerDomTable = HMMER.out[1].collect(), koList = params.koList)
    }

    //CHECKM(DASTOOL.out.dasBins)
    //GTDBTK(DASTOOL.out.dasBins, params.gtdbtkDB)
    //COLLECT(SEQKIT.out, CHECKM.out, GTDBTK.out)

}
