// Mandatory parameters have no default (null)
params.sampleSheet = null
params.odbLocation = null
params.lineage     = null
params.help        = false

def helpMessage() {
    log.info """
    Usage:
      nextflow run main.nf --samplesheet samples.csv --lineage aves --odbLocation /path/to/odb

    Mandatory arguments:
      --sampleSheet   CSV with columns: sample_id,fastqs (multiple FASTQs separated by ';')
      --lineage       BUSCO lineage of organism being assembled (e.g., "aves")
      --odbLocation   location of ODB so it doesn't need to be downloaded each run

    Optional arguments:
      --help          Show this message and exit
    """.stripIndent()
}

def validateParams() {
    def mandatory = ['sampleSheet', 'lineage', 'odbLocation']

    // "--foo" with no value arrives as boolean true, so treat that as missing too
    def missing = mandatory.findAll { name ->
        def v = params[name]
        v == null || v == true || v.toString().trim().isEmpty()
    }

    if (missing) {
        helpMessage()
        error "Missing mandatory parameter(s): ${missing.collect { '--' + it }.join(', ')}"
    }

    if (!file(params.samplesheet).exists())
        error "Sample sheet not found: ${params.samplesheet}"
    if (!file(params.odbLocation).exists())
        error "ODB not found: ${params.reference}"
}

process hifiasm {
    conda 'hifiasm=0.25.0'

    input:
        tuple val(sample_id), path(fastqs)

    output:
        tuple val(sample_id), path("*.p_ctg.gfa", arity: '1..*')

    script:
    """
    hifiasm -t ${task.cpus} -o ${sample_id} ${fastqs}
    """
}

process gfa2fa {
    input:
        tuple val(sample_id), path(gfa)

    output:
        tuple val(sample_id), path("*.fa")

    script:
    """
    awk '\$1=="S" {print ">" \$2 "\\n" \$3}' $gfa > ${gfa.baseName}.fa
    """
}

process assemblyStats {
    conda 'assembly-stats=1.0.1'

    input:
        tuple val(sample_id), path(fa)

    output:
        tuple val(sample_id), path("*.stats")

    script:
    """
    assembly-stats -s $fa > ${fa.baseName}.stats
    """
}

process compleasm {
    conda 'compleasm=0.2.10'

    input:
        tuple val(sample_id), path(fa)

    output:
        tuple val(sample_id), path("*-compleasm/")

    script:
    """
    compleasm run -t ${task.cpus} \
        -l ${params.lineage} \
        -L ${params.odbLocation} \
        -a $fa \
        -o ${fa.baseName}-compleasm
    """
}

process inspector {
    conda 'inspector=1.3.1'

    input:
        tuple val(sample_id), path(contigs), path(reads)

    output:
        tuple val(sample_id), path("*.inspector_stats")

    script:
    """
    inspector.py -t ${task.cpus} \
        -d pacbio-hifi \
        -r $reads \
        -c $contigs \
        -o inspector_out

    mv inspector_out/summary_statistics ${contigs.baseName}.inspector_stats
    """
}

workflow {
    main:
    if (params.help) {
        helpMessage()
        exit 0
    }
    validateParams()

    samples = Channel
        .from(params.sampleSheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, row.fastqs.tokenize(';').collect { file(it) }) }
    hifiasm(samples)
    
    assemblies = hifiasm.out.transpose() | gfa2fa
    assemblyStats(assemblies)
    compleasm(assemblies)
    inspector(assemblies.combine(samples, by: 0))

    publish:
    assemblies_out = assemblies
    assembly_stats = assemblyStats.out
    busco_scores = compleasm.out
    qv_stats = inspector.out
}

output {
    assemblies {
        path { "${it[0]}/assemblies/"}
    }

    assembly_stats {
        path { "${it[0]}/qc/" }
    }

    busco_scores {
        path { "${it[0]}/qc/" }
    }

    qv_stats {
        path { "${it[0]}/qc/" }
    }
}
