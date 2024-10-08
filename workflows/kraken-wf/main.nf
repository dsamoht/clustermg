include { BRACKEN                           } from '../../modules/bracken'
//include { CONCATENATE_FASTQ                 } from '../../modules/concatenate_fastq'
include { KRAKEN                            } from '../../modules/kraken'
include { KRAKENTOOLS_KRONA                 } from '../../modules/krakentools/krakentools_krona'
include { KRONA                             } from '../../modules/krona'


workflow KRAKEN_WF {
    """
    Run Kraken on set of reads determined by the user
    (`--krakenReads` argument). Paired-end reads must
    merged into a single file.
    """
    take:
    long_reads
    short_reads

    main:
    //if (params.krakenReads == ""){
    //    params.krakenReads = params.longReads
    //}
    //newest_reads = Channel.watchPath("$projectDir/test_data/*.f*.gz").view()
    //reads = Channel.fromPath("$projectDir/test_data/*.f*.gz")
    //CONCATENATE_FASTQ(newest_reads, reads)
    if (params.longReads != '') {
        read_type = "long"
        KRAKEN(long_reads, params.krakenDB, read_type)
    } else {
        read_type = "paired"
        KRAKEN(short_reads, params.krakenDB, read_type)
    }
    BRACKEN(KRAKEN.out.krakenOutputFile, params.krakenDB)
    KRAKENTOOLS_KRONA(BRACKEN.out.brackenOutputForKrona)
    KRONA(KRAKENTOOLS_KRONA.out)

}
