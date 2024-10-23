import org.yaml.snakeyaml.Yaml
import java.nio.file.Files

def get_publishMode(d, mode) {
    def req_mode = mode

    if (req_mode != "auto" && req_mode != "link") {
        return mode
    }

    // default to copy
    mode = "copy"

    file(d).mkdirs()

    testFile = file(workflow.workDir + "/.test")
    testFile.write("test")

    testLink = file(d + "/.test")

    // let's see if we can create hard links
    try {
        Files.createLink(testLink, testFile);
        mode = "link"
    } catch (IOException e) {
        if (req_mode == "link") {
            System.err.println("WARNING: using copy as publish mode, reason: " + e)
        }
    }

    testLink.delete()
    testFile.delete()

    return mode
}

publishDirMode = get_publishMode(params.outputDir, params.publishDirMode)

padding = params.readLength + 100

// Handle BAM input files. We need to convert BAMs to Fastq

process check_PE {
    label 'nextNEOpiENV'
    tag "$meta.sampleName : $meta.sampleType"

    input:
    tuple val(meta), path(bam)
    val nextNEOpiENV_setup

    output:
    tuple val(meta), path(bam), stdout

    script:
    """
    #! /bin/bash
    source \$(conda info --json | awk '/conda_prefix/ { gsub(/"|,/, "", \$2); print \$2 }')/bin/activate /home/dev/manual_conda_envs/nextNEOpiENV
    check_pe.py $bam
    """
}


process bam2fastq {
    label 'nextNEOpiENV'
    tag "$meta.sampleName : $meta.sampleType"
    publishDir "${params.outputDir}/analyses/${meta.sampleName}/01_preprocessing",
        mode: publishDirMode

    input:
    tuple val(meta), path(bam), val(libType), val(libOK)
    val nextNEOpiENV_setup

    output:
    tuple val(meta), path("${prefix}_R*.fastq.gz")

    script:
    def STperThreadMem = (int) Math.max(((int) Math.floor((task.memory.toGiga() - 4) / task.cpus)), 1)
    prefix = meta.sampleName + "_" + meta.sampleType
    meta.libType = libType
    if (libType == "PE")
        """
        #! /bin/bash
        source \$(conda info --json | awk '/conda_prefix/ { gsub(/"|,/, "", \$2); print \$2 }')/bin/activate /home/dev/manual_conda_envs/nextNEOpiENV

        samtools sort -@ ${task.cpus} -m ${STperThreadMem}G -l 0 -n ${bam} | \\
        samtools fastq \\
            -@ ${task.cpus} \\
            -c 5 \\
            -1 ${prefix}_R1.fastq.gz \\
            -2 ${prefix}_R2.fastq.gz \\
            -0 /dev/null -s /dev/null \\
            -n \\
            /dev/stdin
        """
    else if (libType == "SE")
        """
        samtools fastq \\
            -@ ${task.cpus} \\
            -n \\
            ${bam} | \\
            bgzip -@ ${task.cpus} -c /dev/stdin > ${prefix}_R1.fastq.gz
        """
}

process merge_fastq {

    tag "$meta.sampleName : $meta.sampleType"

    publishDir "${params.outputDir}/${meta.sampleName}/01_preprocessing/",
        mode: publishDirMode,
        enabled: params.fullOutput

    input:
    tuple val(meta), path("R1/R1.*"), path("R2/R2.*")

    output:
    tuple val(meta), path("*_R{1,2}.merged.fastq.gz")

    script:
    def prefix = meta.sampleName + "_" + meta.sampleType
    """
    cat R1/* > ${prefix}_R1.merged.fastq.gz
    cat R2/* > ${prefix}_R2.merged.fastq.gz
    """
}

process regions_bed_to_interval_list {

    label 'nextNEOpiENV'
    tag 'RegionsBedToIntervalList'
    publishDir "$params.outputDir/supplemental/00_prepare_Intervals/",
        mode: publishDirMode

    input:
    tuple path(RefDict), path(RegionsBed)
    val nextNEOpiENV_setup

    output:
    path("${RegionsBed.baseName}.interval_list")
    script:
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    """
    #! /bin/bash
    source \$(conda info --json | awk '/conda_prefix/ { gsub(/"|,/, "", \$2); print \$2 }')/bin/activate /home/dev/manual_conda_envs/nextNEOpiENV
    gatk --java-options ${JAVA_Xmx} BedToIntervalList \\
        -I ${RegionsBed} \\
        -O ${RegionsBed.baseName}.interval_list \\
        -SD $RefDict
    """
}

process baits_bed_to_interval_list {
    label 'nextNEOpiENV'
    tag 'BaitsBedToIntervalList'
    publishDir "$params.outputDir/supplemental/00_prepare_Intervals/",
        mode: publishDirMode
    input:
    tuple path(RefDict), path(BaitsBed)
    val nextNEOpiENV_setup
    output:
    path("${BaitsBed.baseName}.interval_list")
    script:
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    """
    #! /bin/bash
    source \$(conda info --json | awk '/conda_prefix/ { gsub(/"|,/, "", \$2); print \$2 }')/bin/activate /home/dev/manual_conda_envs/nextNEOpiENV
    gatk --java-options ${JAVA_Xmx} BedToIntervalList \\
        -I ${BaitsBed} \\
        -O ${BaitsBed.baseName}.interval_list \\
        -SD $RefDict
    """
}

process preprocess_interval_list {

    label 'nextNEOpiENV'

    tag 'preprocessIntervalList'

    publishDir "$params.outputDir/supplemental/00_prepare_Intervals/",
        mode: publishDirMode

    input:
    tuple path(RefFasta), path(RefIdx), path(RefDict)
    path(interval_list)
    val nextNEOpiENV_setup

    output:
    path(outFileName)

    script:
    outFileName = (params.WES) ? interval_list.baseName + "_merged_padded.interval_list" : "wgs_ScatterIntervalsByNs.interval_list"
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    if(params.WES){
        """
        #! /bin/bash
        source \$(conda info --json | awk '/conda_prefix/ { gsub(/"|,/, "", \$2); print \$2 }')/bin/activate /home/dev/manual_conda_envs/nextNEOpiENV
        gatk PreprocessIntervals \\
            -R ${RefFasta} \\
            -L ${interval_list} \\
            --bin-length 0 \\
            --padding ${padding} \\
            --interval-merging-rule OVERLAPPING_ONLY \\
            -O ${outFileName}
        """
    } else {
        """
        #! /bin/bash
        source \$(conda info --json | awk '/conda_prefix/ { gsub(/"|,/, "", \$2); print \$2 }')/bin/activate /home/dev/manual_conda_envs/nextNEOpiENV
        gatk --java-options ${JAVA_Xmx} ScatterIntervalsByNs \\
            --REFERENCE $RefFasta \\
            --OUTPUT_TYPE ACGT \\
            --OUTPUT ${outFileName}
        """
    }
}

process split_intervals {

    label 'nextNEOpiENV'

    tag "SplitIntervals"

    publishDir "$params.outputDir/supplemental/00_prepare_Intervals/SplitIntervals/",
        mode: publishDirMode,
        enabled: params.fullOutput

    input:
    tuple path(RefFasta), path(RefIdx), path(RefDict)
    path(IntervalsList)
    val x
    val nextNEOpiENV_setup

    output:
    path("${IntervalName}/*-scattered.interval_list"), emit: path_list
    val("${IntervalName}"), emit: interval_name

    script:
    IntervalName = IntervalsList.baseName
    """
    #! /bin/bash
    source \$(conda info --json | awk '/conda_prefix/ { gsub(/"|,/, "", \$2); print \$2 }')/bin/activate /home/dev/manual_conda_envs/nextNEOpiENV
    mkdir -p ${params.tmpDir}

    gatk SplitIntervals \\
        --tmp-dir ${params.tmpDir} \\
        -R ${RefFasta}  \\
        -scatter ${x} \\
        --interval-merging-rule ALL \\
        --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \\
        -L ${IntervalsList} \\
        -O ${IntervalName}

    """
}

process interval_list_to_bed {

    label 'nextNEOpiENV'

    tag 'BedFromIntervalList'

    publishDir "$params.outputDir/supplemental/00_prepare_Intervals/",
        mode: publishDirMode

    input:
    path(paddedIntervalList)
    val nextNEOpiENV_setup

    output:
    path("${paddedIntervalList.baseName}.{bed.gz,bed.gz.tbi}")

    script:
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    """
    #! /bin/bash
    source \$(conda info --json | awk '/conda_prefix/ { gsub(/"|,/, "", \$2); print \$2 }')/bin/activate /home/dev/manual_conda_envs/nextNEOpiENV
    gatk --java-options ${JAVA_Xmx} IntervalListToBed \\
        -I ${paddedIntervalList} \\
        -O ${paddedIntervalList.baseName}.bed

    bgzip -c ${paddedIntervalList.baseName}.bed > ${paddedIntervalList.baseName}.bed.gz &&
    tabix -p bed ${paddedIntervalList.baseName}.bed.gz
    """
}

process scattered_interval_list_to_bed {

    label 'nextNEOpiENV'

    tag 'ScatteredIntervalListToBed'

    publishDir "$params.outputDir/supplemental/00_prepare_Intervals/SplitIntervals/${IntervalName}",
        mode: publishDirMode,
        enabled: params.fullOutput

    input:
    tuple val(IntervalName), file(IntervalsList)
    val nextNEOpiENV_setup


    output:
    path("*.bed")

    script:
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    """
    #! /bin/bash
    source \$(conda info --json | awk '/conda_prefix/ { gsub(/"|,/, "", \$2); print \$2 }')/bin/activate /home/dev/manual_conda_envs/nextNEOpiENV
    gatk --java-options ${JAVA_Xmx} IntervalListToBed \\
        -I ${IntervalsList} \\
        -O ${IntervalsList.baseName}.bed
    """
}

process fast_qc {

    label 'fastqc'

    tag "$meta.sampleName : $meta.sampleType"

    publishDir "${params.outputDir}/analyses/${meta.sampleName}/QC/fastqc",
        mode: publishDirMode,
        saveAs: { filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}


    input:
    tuple val(meta), path(reads)


    output:
    tuple val(meta), path("*_fastqc.zip")

    script:
    def reads_R1_ext = (reads[0].getExtension() == "gz") ? "fastq.gz" : reads[0].getExtension()
    def reads_R1     = meta.sampleName + "_" + meta.sampleType + "_R1." + reads_R1_ext

    // do we have PE reads?
    def reads_R2 = "_missing_"
    if(meta.libType == "PE") {
        def reads_R2_ext = (reads[1].getExtension() == "gz") ? "fastq.gz" : reads[1].getExtension()
        reads_R2     = meta.sampleName + "_" + meta.sampleType + "_R2." + reads_R2_ext
    }
    """
    if [ ! -e ${reads_R1} ]; then
        ln -s ${reads[0]} ${reads_R1}
    fi

    if [ "${reads_R2}" != "_missing_" ] && [ ! -e ${reads_R2} ]; then
        ln -s ${reads[1]} ${reads_R2}
    fi

    fastqc --quiet --threads ${task.cpus} \\
        ${reads_R1} ${reads_R2}
    """
}

process fastp {
    label 'fastp'
    tag "$meta.sampleName : $meta.sampleType"
    publishDir "$params.outputDir/analyses/${meta.sampleName}/",
        mode: publishDirMode,
        saveAs: {
            filename ->
                if(filename.indexOf(".json") > 0) {
                    return "QC/fastp/$filename"
                } else if(filename.indexOf("NO_FILE") >= 0) {
                    return null
                } else {
                    return  "01_preprocessing/$filename"
                }
        }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed_R{1,2}.fastq.gz"), emit: reads
    tuple val(meta), path("*.json"), emit: json

    script:
    def reads_R1         = "--in1 " + reads[0]
    def trimmed_reads_R1 = "--out1 " + meta.sampleName + "_" + meta.sampleType + "_trimmed_R1.fastq.gz"
    // do we have PE reads?
    def reads_R2         = ""
    def trimmed_reads_R2 = ""
    if(meta.libType == "PE") {
        reads_R2          = "--in2 " + reads[1]
        trimmed_reads_R2  = "--out2 " + meta.sampleName + "_" + meta.sampleType + "_trimmed_R2.fastq.gz"
    }
    def fastpAdapter = ''
    def adapterSeqFile
    def aseq = false
    def aseqR2 = false
    def afile = false
    if (meta.sampleType.indexOf("DNA") > 0) {
        afile = params.adapterSeqFile
        aseq = params.adapterSeq
        aseqR2 = params.adapterSeqR2
    } else {
         afile = params.adapterSeqFileRNAseq
         aseq = params.adapterSeqRNAseq
         aseqR2 = params.adapterSeqR2RNAseq
    }
    if(afile != false) {
        adapterSeqFile = Channel.fromPath(afile)
        fastpAdapter = "--adapter_fasta " + adapterSeqFile
    } else {
        if(aseq != false) {
            adapterSeq   = Channel.value(aseq)
            fastpAdapter = "--adapter_sequence " + aseq.getVal()
            if(aseqR2 != false && meta.libType == "PE") {
                adapterSeqR2   = Channel.value(aseqR2)
                fastpAdapter += " --adapter_sequence_r2 " + adapterSeqR2.getVal()
            }
        }
    }
    """
    fastp --thread ${task.cpus} \\
        ${reads_R1} \\
        ${reads_R2} \\
        ${trimmed_reads_R1} \\
        ${trimmed_reads_R2} \\
        --json ${meta.sampleName}_${meta.sampleType}_fastp.json \\
        ${fastpAdapter} \\
        ${params.fastpOpts}
    """
}


process fast_qc_trimmed {
    label 'fastqc'
    tag "$meta.sampleName : $meta.sampleType"
    publishDir "${params.outputDir}/analyses/${meta.sampleName}/QC/fastqc",
        mode: publishDirMode,
        saveAs: { filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_fastqc.zip")

    script:
    def reads_R1 = reads[0]
    def reads_R2 = (meta.libType == "PE") ? reads[1] : ""
    """
    fastqc --quiet --threads ${task.cpus} \\
        ${reads_R1} ${reads_R2}
    """
}

process make_ubam {

    label 'nextNEOpiENV'

    tag "$meta.sampleName : $meta.sampleType"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/01_preprocessing/",
        mode: publishDirMode

    input:
    tuple val(meta), path(reads)

    val nextNEOpiENV_setup

    output:
    tuple val(meta), path(ubam)

    script:
    ubam = meta.sampleName + "_" + meta.sampleType + "_unaligned.bam"
    def read_group = meta.sampleName + "_" + meta.sampleType.replaceAll("_DNA", "")
    def reads_in = "-F1 " + reads[0]
    reads_in += (meta.libType == "PE") ? " -F2 " + reads[1] : ""
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    def java_opts = '"' + JAVA_Xmx + ' -XX:ParallelGCThreads=' + task.cpus + '"'

    """
    #! /bin/bash
    source \$(conda info --json | awk '/conda_prefix/ { gsub(/"|,/, "", \$2); print \$2 }')/bin/activate /home/dev/manual_conda_envs/nextNEOpiENV
    mkdir -p ${params.tmpDir}
    gatk --java-options ${java_opts} FastqToSam \\
        --TMP_DIR ${params.tmpDir} \\
        --MAX_RECORDS_IN_RAM ${params.maxRecordsInRam} \\
        ${reads_in} \\
        --READ_GROUP_NAME ${read_group} \\
        --SAMPLE_NAME ${read_group} \\
        --LIBRARY_NAME ${read_group} \\
        --PLATFORM ILLUMINA \\
        -O ${ubam}
    """
}

process bwa {

    label 'nextNEOpiENV'

    tag "$meta.sampleName : $meta.sampleType"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/02_alignments/",
        mode: publishDirMode

    input:
    tuple val(meta), path(reads)

    tuple path(RefFasta), path(RefIdx), path(RefDict), path(BwaRef)

    val nextNEOpiENV_setup

    output:
    tuple val(meta), path(bam)

    script:
    bam = meta.sampleName + "_" + meta.sampleType + "_aligned.bam"
    def read_group = meta.sampleName + "_" + meta.sampleType.replaceAll("_DNA", "")

    def sort_threads = (task.cpus.compareTo(8) == 1) ? 8 : task.cpus
    def SB_sort_mem =  Math.max((task.memory.toGiga() - 4), 1) + "G"
    """
    #! /bin/bash
    source \$(conda info --json | awk '/conda_prefix/ { gsub(/"|,/, "", \$2); print \$2 }')/bin/activate /home/dev/manual_conda_envs/nextNEOpiENV
    bwa mem \\
        -R "@RG\\tID:${read_group}\\tLB:${read_group}\\tSM:${read_group}\\tPL:ILLUMINA" \\
        -M ${RefFasta} \\
        -t ${task.cpus} \\
        -Y \\
        ${reads.join(" ")} | \\
    samtools view -@2 -Shbu - | \\
    sambamba sort \\
        --sort-picard \\
        --tmpdir=${params.tmpDir} \\
        -m ${SB_sort_mem} \\
        -l 6 \\
        -t ${sort_threads} \\
        -o ${bam} \\
        /dev/stdin
    """
}

process merge_ubam_bam {

    label 'nextNEOpiENV'

    tag "$meta.sampleName : $meta.sampleType"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/02_alignments/",
        mode: publishDirMode,
        enabled: params.fullOutput

    input:
    tuple val(meta), path(bam), path(ubam)

    tuple path(RefFasta), path(RefIdx), path(RefDict)

    val nextNEOpiENV_setup

    output:
    tuple val(meta), path("${procSampleName}_aligned_uBAM_merged.bam")

    script:
    procSampleName = meta.sampleName + "_" + meta.sampleType

    def paired_run = (meta.libType == 'SE') ? 'false' : 'true'
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    def java_opts = '"' + JAVA_Xmx + ' -XX:ParallelGCThreads=' + task.cpus + '"'
    """
    #! /bin/bash
    source \$(conda info --json | awk '/conda_prefix/ { gsub(/"|,/, "", \$2); print \$2 }')/bin/activate /home/dev/manual_conda_envs/nextNEOpiENV
    mkdir -p ${params.tmpDir}

    gatk --java-options ${java_opts} MergeBamAlignment \\
        --TMP_DIR ${params.tmpDir} \\
        --VALIDATION_STRINGENCY SILENT \\
        --EXPECTED_ORIENTATIONS FR \\
        --ATTRIBUTES_TO_RETAIN X0 \\
        --REFERENCE_SEQUENCE ${RefFasta} \\
        --PAIRED_RUN ${paired_run} \\
        --SORT_ORDER "queryname" \\
        --IS_BISULFITE_SEQUENCE false \\
        --ALIGNED_READS_ONLY false \\
        --CLIP_ADAPTERS false \\
        --MAX_RECORDS_IN_RAM ${params.maxRecordsInRamMerge} \\
        --ADD_MATE_CIGAR true \\
        --MAX_INSERTIONS_OR_DELETIONS -1 \\
        --PRIMARY_ALIGNMENT_STRATEGY MostDistant \\
        --UNMAPPED_READ_STRATEGY COPY_TO_TAG \\
        --ALIGNER_PROPER_PAIR_FLAGS true \\
        --UNMAP_CONTAMINANT_READS true \\
        --ALIGNED_BAM ${bam} \\
        --UNMAPPED_BAM ${ubam} \\
        --OUTPUT ${procSampleName}_aligned_uBAM_merged.bam
    """
}

process mark_duplicates {

    label 'nextNEOpiENV'

    tag "$meta.sampleName : $meta.sampleType"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/02_alignments/",
        mode: publishDirMode

    input:
    tuple val(meta), path(bam)

    tuple path(RefFasta), path(RefIdx), path(RefDict)

    val nextNEOpiENV_setup

    output:
    tuple val(meta), path(bam_out)
    script:
    def procSampleName = meta.sampleName + "_" + meta.sampleType
    def STperThreadMem = (int) Math.max(((int) Math.floor((task.memory.toGiga() - 8) / task.cpus)), 1)
    def JAVA_Xmx = '-Xmx4G'
    bam_out = [procSampleName + "_aligned_sort_mkdp.bam", procSampleName + "_aligned_sort_mkdp.bai"]
    """
    #! /bin/bash
    source \$(conda info --json | awk '/conda_prefix/ { gsub(/"|,/, "", \$2); print \$2 }')/bin/activate /home/dev/manual_conda_envs/nextNEOpiENV
    mkdir -p ${params.tmpDir}
    sambamba markdup \\
        -t ${task.cpus} \\
        --tmpdir ${params.tmpDir} \\
        --hash-table-size=${params.SB_hash_table_size } \\
        --overflow-list-size=${params.SB_overflow_list_size} \\
        --io-buffer-size=${params.SB_io_buffer_size} \\
        ${bam} \\
        /dev/stdout | \\
    samtools sort \\
        -@${task.cpus} \\
        -m ${STperThreadMem}G \\
        -O BAM \\
        -l 0 \\
        /dev/stdin | \\
    gatk --java-options ${JAVA_Xmx} SetNmMdAndUqTags \\
        --TMP_DIR ${params.tmpDir} \\
        -R ${RefFasta} \\
        -I /dev/stdin \\
        -O ${procSampleName}_aligned_sort_mkdp.bam \\
        --CREATE_INDEX true \\
        --MAX_RECORDS_IN_RAM ${params.maxRecordsInRam} \\
        --VALIDATION_STRINGENCY LENIENT

    """
}

process alignment_metrics {

    label 'nextNEOpiENV'
    tag "$meta.sampleName : $meta.sampleType"
    publishDir "$params.outputDir/analyses/${meta.sampleName}/QC/alignments/",
        mode: publishDirMode

    input:
    tuple val(meta), path(bam)
    tuple path(RefFasta), path(RefIdx)
    path(BaitIntervalsList)
    path(IntervalsList)
    val nextNEOpiENV_setup

    output:
    tuple val(meta), path("${procSampleName}.*.txt")

    script:
    procSampleName = meta.sampleName + "_" + meta.sampleType
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    def java_opts = '"' + JAVA_Xmx + ' -XX:ParallelGCThreads=' + task.cpus + '"'
    """
    #! /bin/bash
    source \$(conda info --json | awk '/conda_prefix/ { gsub(/"|,/, "", \$2); print \$2 }')/bin/activate /home/dev/manual_conda_envs/nextNEOpiENV
    mkdir -p ${params.tmpDir}
    gatk --java-options ${java_opts} CollectHsMetrics \\
        --TMP_DIR ${params.tmpDir} \\
        --INPUT ${bam[0]} \\
        --OUTPUT ${procSampleName}.HS.metrics.txt \\
        -R ${RefFasta} \\
        --BAIT_INTERVALS ${BaitIntervalsList} \\
        --TARGET_INTERVALS ${IntervalsList} \\
        --PER_TARGET_COVERAGE ${procSampleName}.perTarget.coverage.txt && \\
    gatk --java-options ${java_opts} CollectAlignmentSummaryMetrics \\
        --TMP_DIR ${params.tmpDir} \\
        --INPUT ${bam[0]} \\
        --OUTPUT ${procSampleName}.AS.metrics.txt \\
        -R ${RefFasta} &&
    samtools flagstat -@${task.cpus} ${bam[0]} > ${procSampleName}.flagstat.txt
    """
}

process scatter_base_recal_gatk4 {

    label 'nextNEOpiENV'

    tag "$meta.sampleName : $meta.sampleType"

    input:
    tuple val(meta), path(bam), path(intervals)

    tuple path(RefFasta), path(RefIdx), path(RefDict)

    tuple path(MillsGold), path(MillsGoldIdx), path(DBSNP), path(DBSNPIdx), path(KnownIndels), path(KnownIndelsIdx)

    val nextNEOpiENV_setup

    output:
    tuple val(meta), path("${procSampleName}_${intervals}_bqsr.table")


    script:
    procSampleName = meta.sampleName + "_" + meta.sampleType
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    """
    #! /bin/bash
    source \$(conda info --json | awk '/conda_prefix/ { gsub(/"|,/, "", \$2); print \$2 }')/bin/activate /home/dev/manual_conda_envs/nextNEOpiENV
    mkdir -p ${params.tmpDir}
    gatk  --java-options ${JAVA_Xmx} BaseRecalibrator \\
        --tmp-dir ${params.tmpDir} \\
        -I ${bam[0]} \\
        -R ${RefFasta} \\
        -L ${intervals} \\
        -O ${procSampleName}_${intervals}_bqsr.table \\
        --known-sites ${DBSNP} \\
        --known-sites ${KnownIndels} \\
        --known-sites ${MillsGold}
    """
}

process gather_gatk4_scsattered_bqsr_tables {

    label 'nextNEOpiENV'

    tag "$meta.sampleName : $meta.sampleType"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/03_baserecalibration/",
        mode: publishDirMode,
        enabled: params.fullOutput

    input:
    tuple val(meta), path(bqsr_table)
    val nextNEOpiENV_setup

    output:
    tuple val(meta), path("${procSampleName}_bqsr.table")

    script:
    procSampleName = meta.sampleName + "_" + meta.sampleType

    """
    #! /bin/bash
    source \$(conda info --json | awk '/conda_prefix/ { gsub(/"|,/, "", \$2); print \$2 }')/bin/activate /home/dev/manual_conda_envs/nextNEOpiENV
    mkdir -p ${params.tmpDir}

    gatk GatherBQSRReports \\
        -I ${bqsr_table.join(" -I ")} \\
        -O ${procSampleName}_bqsr.table
    """
}

process scatter_gatk4_apply_bqsrs {

    label 'nextNEOpiENV'

    tag "$meta.sampleName : $meta.sampleType"

    input:
    tuple val(meta), path(bam), path(bqsr_table), path(intervals)

    tuple path(RefFasta), path(RefIdx), path(RefDict)

    tuple path(MillsGold), path(MillsGoldIdx), path(DBSNP), path(DBSNPIdx), path(KnownIndels), path(KnownIndelsIdx)

    val nextNEOpiENV_setup

    output:
    tuple val(meta), path(bam_out)


    script:
    def procSampleName = meta.sampleName + "_" + meta.sampleType
    bam_out = [ procSampleName + "_" + intervals + "_recal4.bam",
                procSampleName + "_" + intervals + "_recal4.bai"]
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    """
    #! /bin/bash
    source \$(conda info --json | awk '/conda_prefix/ { gsub(/"|,/, "", \$2); print \$2 }')/bin/activate /home/dev/manual_conda_envs/nextNEOpiENV
    mkdir -p ${params.tmpDir}
    gatk ApplyBQSR \\
        --java-options ${JAVA_Xmx} \\
        --tmp-dir ${params.tmpDir} \\
        -I ${bam[0]} \\
        -R ${RefFasta} \\
        -L ${intervals} \\
        -O ${procSampleName}_${intervals}_recal4.bam \\
        --bqsr-recal-file ${bqsr_table}
    """
}

process gather_recal_bam_files {

    label 'nextNEOpiENV'

    tag "$meta.sampleName : $meta.sampleType"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/03_baserecalibration/",
        mode: publishDirMode

    input:
    tuple val(meta), path(bam), path(bai)

    val nextNEOpiENV_setup

    output:
    tuple val(meta), path("${procSampleName}_recalibrated.{bam,bam.bai}")

    script:
    procSampleName = meta.sampleName + "_" + meta.sampleType
    def STperThreadMem = (int) Math.max(((int) Math.floor((task.memory.toGiga() - 4) / task.cpus)), 1)
    def JAVA_Xmx = "4G"
    def java_opts = '"-Xmx' + JAVA_Xmx + ' -XX:ParallelGCThreads=2"'
    """
    #! /bin/bash
    source \$(conda info --json | awk '/conda_prefix/ { gsub(/"|,/, "", \$2); print \$2 }')/bin/activate /home/dev/manual_conda_envs/nextNEOpiENV
    mkdir -p ${params.tmpDir}

    rm -f ${procSampleName}_gather.fifo
    mkfifo ${procSampleName}_gather.fifo
    gatk --java-options ${java_opts} GatherBamFiles \\
        --TMP_DIR ${params.tmpDir} \\
        -I ${bam.join(" -I ")} \\
        -O ${procSampleName}_gather.fifo \\
        --CREATE_INDEX false \\
        --MAX_RECORDS_IN_RAM ${params.maxRecordsInRam} &
    samtools sort \\
        -@${task.cpus} \\
        -m ${STperThreadMem}G \\
        -o ${procSampleName}_recalibrated.bam ${procSampleName}_gather.fifo
    samtools index -@${task.cpus} ${procSampleName}_recalibrated.bam
    rm -f ${procSampleName}_gather.fifo
    """
}
