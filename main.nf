nextflow.enable.dsl=2

import org.yaml.snakeyaml.Yaml
import java.nio.file.Files

include {
  check_PE; bam2fastq; merge_fastq; regions_bed_to_interval_list;
  baits_bed_to_interval_list; preprocess_interval_list; split_intervals;
  interval_list_to_bed; scattered_interval_list_to_bed; fast_qc;
  fastp; fast_qc_trimmed; make_ubam; bwa; merge_ubam_bam; mark_duplicates;
  alignment_metrics; scatter_base_recal_gatk4; gather_gatk4_scsattered_bqsr_tables;
  scatter_gatk4_apply_bqsrs; gather_recal_bam_files;
} from './modules/pre-processing.nf'


def checkDir(d, description) {
    myDir = file(d)
    result = myDir.mkdirs()
    if (result) {
        println(description + ": " + myDir.toRealPath())
    } else {
        exit 1, "Cannot create directory: " + myDir
    }
    return myDir.toRealPath()
}

def mkTmpDir(d) {
   return checkDir(d, "tmpDir")
}

def check_iedb_dir(d) {
    checkDir(d, "IEDB_dir")
}

def check_mhcflurry_dir(d) {
    checkDir(d, "MHCFLURRY_dir")
}

def setExomeCaptureKit(captureKit) {
    resources = ['BaitsBed', 'RegionsBed']
    for (r in resources) {
        if (params.exomeCaptureKits[ captureKit ][r] == null) {
            exit 1, "ERROR: \"" + r + "\" file not set for: " + captureKit + "\nPlease check the \"exomeCaptureKits\" resource file settings in conf/resources.config"
        }
    }

    params.references.BaitsBed = params.exomeCaptureKits[ captureKit ].BaitsBed
    params.references.RegionsBed = params.exomeCaptureKits[ captureKit ].RegionsBed
}

def check_resource(resource, resource_type) {

    resource_file = params[resource_type][resource]

    if (resource_file == null) {
        exit 1, "ERROR: Resource file not set for: " + resource + "\nPlease check the \"" + resource_type + "\" resource file settings in conf/resources.config"
    }

    rp = file(resource_file)

    if (rp instanceof java.util.LinkedList) {
        rf = rp
    } else {
        rf = [rp]
    }

    err = 0
    if (rf[0] != null) {
        for (f in rf) {
            err += file(f).exists() ? 0 : 1
        }
    }

    if (err == 0) {
        return(rp)
    } else {
        exit 1, "ERROR: Resource file does not exist: " + resource_file + "\nPlease check the " + resource_type + " resource file settings in conf/resources.config"
    }
}

def defineResources(resource_type, wes, hlahd) {

    // vep check file
    vep_cache_chck_file_name = "." + params.vep_species + "_" + params.vep_assembly + "_" + params.vep_cache_version + "_cache_ok.chck"
    vep_cache_chck_file = file(params.databases.vep_cache + "/" + vep_cache_chck_file_name)

    // define references
    references = ['RefFasta', 'RefIdx', 'RefDict', 'RefChrLen', 'RefChrDir', 'BwaRef',
                  'YaraIndexDNA', 'YaraIndexRNA', 'STARidx', 'AnnoFile', 'ExonsBED', 'acLoci',
                  'acLociGC', 'SequenzaGC', 'ProteinBlastDBdir']
    if (wes) {
        references.addAll(['BaitsBed', 'RegionsBed'])
    }
    if (hlahd != "") {
        references.addAll(['HLAHDFreqData', 'HLAHDGeneSplit', 'HLAHDDict'])
    }
    if(vep_cache_chck_file.exists() && !vep_cache_chck_file.isEmpty()) {
        references.addAll(['VepFasta'])
    }

    // define databases
    databases = ['MillsGold', 'MillsGoldIdx', 'hcSNPS1000G', 'hcSNPS1000GIdx', 'HapMap',
                 'HapMapIdx', 'DBSNP', 'DBSNPIdx', 'GnomAD', 'GnomADIdx', 'GnomADfull', 'GnomADfullIdx',
                 'KnownIndels', 'KnownIndelsIdx', 'vep_cache', 'IEDB_dir', 'MHCFLURRY_dir']

    resources = ['references' : references, 'databases' : databases ]
    resources_files = [:]

    for (r in resources[resource_type]) {
        resources_files[r] = check_resource(r, resource_type)
    }

    return(resources_files)
}

def checkToolAvailable(tool, check, errMode, module=false) {
    def checkResult = false
    def res = ""
    def chckCmd = ""
    def envlist = [];

    if (check == "inPath") {
        if(module) {
            chckCmd = 'module load ' + module + ' && which ' + tool + ' && module unload ' + module
        } else {
            chckCmd = "which " + tool
        }
        def processBuilder = new ProcessBuilder(['/bin/bash', '-c', chckCmd])
        processBuilder.environment().putAll(System.getenv())

        def process = processBuilder.start()
        def out = new StringBuilder()

        process.inputStream.eachLine { line ->
            out.append(line).append('\n')
        }

        process.waitFor()
        res = out.toString().trim()
    }

    if (check == "exists") {
        if (file(tool).exists()) {
            res = tool
        }
    }

    if (res == "") {
        def msg = tool + " not found, please make sure " + tool + " is installed"

        if(errMode == "err") {
            msg = "ERROR: " + msg
            msg = (check == "inPath") ? msg + " and in your \$PATH" : msg
            exit(1, msg)
        } else {
            msg = "Warning: " + msg
            msg = (check == "inPath") ? msg + " and in your \$PATH" : msg
            println("Warning: " + msg)
        }
    } else {
        println("Found " + tool + " at: " + res)
        checkResult = true
    }

    return checkResult
}

def checkCondaChannels() {
    Yaml parser = new Yaml()
    def channels = []
    try {
        def config = parser.load("conda config --show channels".execute().text)
        channels = config.channels
    } catch(NullPointerException | IOException e) {
        log.warn "Could not verify conda channel configuration."
        return
    }

    // Check that all channels are present
    def required_channels = ['conda-forge', 'bioconda', 'defaults']
    def conda_check_failed = !required_channels.every { ch -> ch in channels }

    // Check that they are in the right order
    conda_check_failed |= !(channels.indexOf('conda-forge') < channels.indexOf('bioconda'))
    conda_check_failed |= !(channels.indexOf('bioconda') < channels.indexOf('defaults'))

    if (conda_check_failed) {
        log.warn "=============================================================================\n" +
            "  There is a problem with your Conda configuration!\n\n" +
            "  You will need to set-up the conda-forge and bioconda channels correctly.\n" +
            "  Please refer to https://bioconda.github.io/#usage\n" +
            "  NB: The order of the channels matters!\n" +
            "==================================================================================="

        exit 1
    }
}

def showLicense() {

    licenseFile = file(baseDir + "/LICENSE")
    log.info licenseFile.text

    log.info ""
    log.warn "To accept the licence terms, please rerun with '--accept_license'"
    log.info ""

    exit 1
}

def acceptLicense() {
    log.info ""
    log.warn "I have read and accept the licence terms"
    log.info ""

    licenseChckFile = file(baseDir + "/.license_accepted.chck")
    licenseChckFile.text = "License accepted by " + workflow.userName + " on "  + workflow.start

    return true
}

def checkLicense() {
    licenseChckFile = file(baseDir + "/.license_accepted.chck")

    if(!licenseChckFile.exists()) {
        showLicense()
    } else {
        return true
    }
}


def check_seqLibTypes_ok(seqLib_ch) {
    def seqLibs = seqLib_ch.toList().get()
    def pe_count = 0
    def se_count = 0

    def lt_map = [:]

    for (seqLib in seqLibs) {
        if (seqLib[0].sampleType != "tumor_RNA") {
            if (! lt_map.containsKey(seqLib[0].sampleType)) {
                lt_map[seqLib[0].sampleName] = seqLib[2]
            } else {
                if (lt_map[seqLib[0].sampleName] != seqLib[2]) {
                    exit 1, "Please do not mix pe and se for tumor/normal pairs: " + seqLib[0].sampleName + " - Not supported"
                }
            }
        }
    }
    return "OK"
}

// This function removes all keys in key list from meta object at idx 0.
// If keep is true then the original meta object is kept at idx 1
// of the channel values. Note with keep true all values [1..-1] will be
// right shifted by 1
/* NOT working or now: need to check
def remove_from_meta(ch, keys=[], keep=false) {
    ch = ch.map {
            meta, f ->
            def meta_new = meta.clone()
            keys.each{ k ->
                meta_new.remove(k)
            }
            if(keep) {
                return [meta_new, meta.clone(), f]
            } else {
                return [meta_new, f]
            }
        }
    return ch
}
*/



def helpMessage() {
    log.info ""
    log.info "----------------------------"
    log.info "--        U S A G E       "
    log.info "----------------------------"
    log.info ""
    log.info ' nextflow run nextNEOpi.nf -config conf/params.config --batchFile <batchfile.csv> -profile [conda|singularity],[cluster] [-resume]'
    log.info ""
    log.info "-----------------------------------------------------------------------------------------------------------------------------------------"
    log.info ""
    log.info ""
    log.info " Mandatory arguments:"
    log.info " --------------------"
    log.info "--batchFile"
    log.info ""
    log.info "CSV-file, T/N reads, and optionally RNAseq reads:"

    log.info "sampleName,reads1,reads2,sampleType,HLAfile,sex"
    log.info "sample1,reads_s1_t_1.fastq.gz,reads_s1_t_2.fastq.gz,tumor_DNA,,female"
    log.info "sample1,reads_s1_n_1.fastq.gz,reads_s1_n_2.fastq.gz,normal_DNA,,female"
    log.info "sample1,reads_s1_r_1.fastq.gz,reads_s1_r_2.fastq.gz,tumor_RNA,,female"
    log.info "sample2,reads_s2_t_1.fastq.gz,reads_s2_t_2.fastq.gz,tumor_DNA,/data/sample2_hla.txt,male"
    log.info "sample2,reads_s2_n_1.fastq.gz,reads_s2_n_2.fastq.gz,normal_DNA,,male"
    log.info "sample2,reads_s2_r_1.fastq.gz,,tumor_RNA,,male"

    log.info "Note: You can not use samples that have mixed single-end and paired-end DNA reads in tumor and normal."
    log.info "Both, tumor and normal DNA library types need to be either SE or PE for a given sample."

    log.info "Note: in the HLAfile coulumn a user suppiled HLA types file may be specified for a given sample"

    log.info "Note: sex can be XX, female or Female, XY, male or Male. If not specified or \"NA\" the gender is inferred from the data"
    log.info ""
    log.info ""
    log.info ""
    log.info ""
    log.info " All references, databases, software should be edited in the resources.config file"
    log.info ""
    log.info " For further options e.g. BAM input see the README.md, the params.config and process.config files"
    log.info "-----------------------------------------------------------------------------------------------------------------------------------------"
}



workflow {
  log.info ""
  log.info " NEXTFLOW ~  version ${workflow.nextflow.version} ${workflow.nextflow.build}"
  log.info "-------------------------------------------------------------------------"
  log.info "          Nextfow NeoAntigen Prediction Pipeline - nextNEOpi    "
  log.info "-------------------------------------------------------------------------"
  log.info ""
  log.info " Features: "
  log.info " - somatic variants from tumor + matched normal samples"
  log.info " - CNV analysis"
  log.info " - tumor muational burden"
  log.info " - class I and class II HLA typing"
  log.info " - gene fusion peptide prediction using RNAseq data"
  log.info " - peptide MHC binding perdiction"
  log.info " - clonality of neoantigens"
  log.info " - expression of neoantigens"
  log.info ""
  log.info "-------------------------------------------------------------------------"
  log.info "C O N F I G U R A T I O N"
  log.info ""
  log.info "Command Line: \t\t " + workflow.commandLine
  log.info "Working Directory: \t " + workflow.workDir
  log.info "Output Directory: \t " + params.outputDir
  log.info ""
  log.info "I N P U T"
  log.info ""
  log.info "batch file: \t\t " + params.batchFile
  log.info ""
  log.info "Please check --help for further instruction"
  log.info "-------------------------------------------------------------------------"

  // Check if License(s) were accepted
  params.accept_license = false

  if (params.accept_license) {
      acceptLicense()
  } else {
      checkLicense()
  }

  /*
  ________________________________________________________________________________

                              C O N F I G U R A T I O N
  ________________________________________________________________________________
  */
  if (params.help) exit 0, helpMessage()

  // switch for enable/disable processes (debug/devel only: use if(params.RUNTHIS) { .... })
  params.RUNTHIS = false

  // default is not to get bams as input data
  bamInput = false

  // initialize RNA tag seq
  have_RNA_tag_seq = params.RNA_tag_seq

  // set and initialize the Exome capture kit
  setExomeCaptureKit(params.exomeCaptureKit)

  // check conda channels
  if (params.enable_conda) {
      checkCondaChannels()
  }

  // check IEDB dir
  check_iedb_dir(params.databases.IEDB_dir)

  // check MHCflurry dir
  check_mhcflurry_dir(params.databases.MHCFLURRY_dir)

  // set and check references and databases
  reference = defineResources('references', params.WES, params.HLAHD_DIR)
  database = defineResources('databases', params.WES, params.HLAHD_DIR)

  // create tmp dir and make sure we have the realpath for it
  tmpDir = mkTmpDir(params.tmpDir)

  /*--------------------------------------------------
    For workflow summary
  ---------------------------------------------------*/
  // Has the run name been specified by the user?
  //  this has the bonus effect of catching both -name and --name
  custom_runName = params.name
  if ( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ) {
      custom_runName = workflow.runName
  }

  // Summary
  def summary = [:]
  summary['Pipeline Name']                 = 'icbi/nextNEOpi'
  summary['Pipeline Version']              = workflow.manifest.version
  summary['Batch file']                    = params.batchFile
  summary['Read length']                   = params.readLength
  summary['Exome capture kit']             = params.exomeCaptureKit
  summary['Fasta Ref']                     = params.references.RefFasta
  summary['MillsGold']                     = params.databases.MillsGold
  summary['hcSNPS1000G']                   = params.databases.hcSNPS1000G
  summary['HapMap']                        = params.databases.HapMap
  summary['Cosmic']                        = params.databases.Cosmic
  summary['DBSNP']                         = params.databases.DBSNP
  summary['GnomAD']                        = params.databases.GnomAD
  summary['GnomADfull']                    = params.databases.GnomADfull
  summary['KnownIndels']                   = params.databases.KnownIndels
  summary['BlastDB']                       = params.references.ProteinBlastDBdir
  summary['priority variant Caller']       = params.primaryCaller
  summary['Mutect 1 and 2 minAD']          = params.minAD
  summary['VarScan min_cov']               = params.min_cov
  summary['VarScan min_cov_tumor']         = params.min_cov_tumor
  summary['VarScan min_cov_normal']        = params.min_cov_normal
  summary['VarScan min_freq_for_hom']      = params.min_freq_for_hom
  summary['VarScan somatic_pvalue']        = params.somatic_pvalue
  summary['VarScan somatic_somaticpvalue'] = params.somatic_somaticpvalue
  summary['VarScan strand_filter']         = params.strand_filter
  summary['VarScan processSomatic_pvalue'] = params.processSomatic_pvalue
  summary['VarScan max_normal_freq']       = params.max_normal_freq
  summary['VarScan min_tumor_freq']        = params.min_tumor_freq
  summary['VarScan min_map_q']             = params.min_map_q
  summary['VarScan min_base_q']            = params.min_base_q
  summary['VEP assembly']                  = params.vep_assembly
  summary['VEP species']                   = params.vep_species
  summary['VEP options']                   = params.vep_options
  summary['Number of scatters']            = params.scatter_count
  summary['Output dir']                    = params.outputDir
  summary['Working dir']                   = workflow.workDir
  summary['TMP dir']                       = tmpDir
  summary['Current home']                  = "$HOME"
  summary['Current user']                  = "$USER"
  summary['Current path']                  = "$PWD"
  summary['Picard maxRecordsInRam']        = params.maxRecordsInRam
  summary['Script dir']                    = workflow.projectDir
  summary['Config Profile']                = workflow.profile


  if(params.email) summary['E-mail Address'] = params.email
  log.info summary.collect { k,v -> "${k.padRight(30)}: $v" }.join("\n")
  log.info "-------------------------------------------------------------------------"

  // End Summary

  // determine the publishDirMode
//   publishDirMode = get_publishMode(params.outputDir, params.publishDirMode)

  // Check if we got a batch file
  params.batchFile = null
  if (params.batchFile == null) {
      log.error "No sample sheet specified, please use --batchFile to pass the sample sheet"
      exit(1)
  }

  def batchCSV = file(params.batchFile).splitCsv(header:true)

  def validFQfields = ["sampleName",
                       "reads1",
                       "reads2",
                       "sampleType",
                      "HLAfile",
                      "sex"]

  def validBAMfields = ["sampleName",
                        "bam",
                        "sampleType",
                        "HLAfile",
                        "sex"]

  def validSampleTypes = ["tumor_DNA", "normal_DNA", "tumor_RNA"]

  if (batchCSV.size() > 0) {

      if (batchCSV[0].keySet().sort() == validFQfields.sort()) {
          bamInput = false
      } else if (batchCSV[0].keySet().sort() == validBAMfields.sort()) {
          bamInput = true
      } else {
          exit 1, "Error: Incorrect fields in batch file, please check your batchFile"
      }
  } else {
      exit 1, "Error: No samples found, please check your batchFile"
  }

  def raw_data = []
  def custom_HLA_data = []

  def t_map = [:]
  def n_map = [:]
  def r_map = [:]
  def lt_map = [:]

  for ( row in batchCSV ) {
      def meta  = [:]
      def reads = []

      meta.sampleName = row.sampleName.toString()
      meta.sampleType = row.sampleType.toString()

      t_map[meta.sampleName] = (t_map[meta.sampleName] == null) ? 0 : t_map[meta.sampleName]
      n_map[meta.sampleName] = (n_map[meta.sampleName] == null) ? 0 : n_map[meta.sampleName]

      if(row.sex) {
          if (row.sex.toLowerCase() in ["xx", "female"]) {
             meta.sex = "XX"
          } else if (row.sex.toLowerCase() in ["xy", "male"]) {
              meta.sex = "XY"
          } else if (row.sex.toLowerCase() in ["none", "na"]) {
              meta.sex = "None"
              println("WARNING: " + row.sampleName + " sex not specified will infer from data")
          } else {
              exit 1, "sex should be one of: XX, xx, XY, xy, Female, female, Male, male, None, none, NA, got: " + row.sex
          }
          meta.maleRef = (meta.sex ==  "XY") ? true : false
      } else {
          println("WARNING: " + row.sampleName + " sex not specified will infer from data")
          meta.sex = "None"
          meta.maleRef = true
      }

      if (row.sampleType in validSampleTypes) {
          t_map[meta.sampleName] += (row.sampleType == "tumor_DNA") ? 1 : 0
          n_map[meta.sampleName] += (row.sampleType == "normal_DNA") ? 1 : 0
      } else {
          exit 1, "Error: Incorrect sampleType [got: " + row.sampleType + " - allowed: " + validSampleTypes + "]: please check your batchFile"
      }

      // remove stuff that is not required in the HLAfile channel
      def meta_min = meta.clone()
      meta_min.keySet().removeAll(['sampleType', 'libType'])
      if (row.HLAfile) {
          custom_HLA_data.add([meta_min, file(row.HLAfile, checkIfExists: true)])
      } else {
          custom_HLA_data.add([meta_min, []])
      }

      if (row.sampleType == "tumor_RNA") {
          meta.have_RNA = true
          r_map[row.sampleName] = true
      }


      if (! row.bam) {
          meta.libType = "SE"
          if (row.reads1) { reads.add(file(row.reads1, checkIfExists: true)) }
          if (row.reads2) {
              reads.add(file(row.reads2, checkIfExists: true))
              meta.libType = "PE"
          }
          if (meta.sampleType != "tumor_RNA") {
              if (lt_map[meta.sampleName] == null) {
                  lt_map[meta.sampleName] = meta.libType
              } else {
                  if (lt_map[meta.sampleName] != meta.libType) {
                      exit 1, "Please do not mix pe and se for tumor/normal pairs: " + meta.sampleName + " - Not supported"
                  }
              }
          }
          raw_data.add([meta, reads])
      } else {
          raw_data.add([meta, file(row.bam, checkIfExists: true)])
      }

  }

  // check if we have T/N DNA pairs for all patients
  raw_data.each {
      record ->
      if (t_map[record[0].sampleName] == 0) {
          exit 1, "NO tumor DNA sample specified for: " + record[0].sampleName
      }
      if (n_map[record[0].sampleName] == 0) {
          exit 1, "NO normal DNA sample specified for: " + record[0].sampleName
      }
  }

  // update meta for RNA
  raw_data.each {
      record ->
      record[0].have_RNA = (r_map[record[0].sampleName]) ? true : false
  }
  custom_HLA_data = custom_HLA_data.each {
      record ->
      record[0].have_RNA = (r_map[record[0].sampleName]) ? true : false
  }.unique()

  use_custom_hlas = (custom_HLA_data.size() > 0) ? true : false

  batch_raw_data_ch = Channel.fromList(raw_data)
  batch_custom_HLA_data_ch = Channel.fromList(custom_HLA_data)

  if (bamInput == false) {
      batch_raw_data_ch.map {
              meta, fastq ->
              [ meta, fastq ]
          }
          .groupTuple(by: [0])
          .branch {
              meta, fastq ->
                  single  : fastq.unique().size() == 1
                      return [ meta, fastq.flatten().unique() ]
                  multiple: fastq.unique().size() > 1
                      return [ meta, fastq.flatten().unique() ]
          }
          .set { fastq_ch }
  } else {
      fastq_ch = Channel.of().branch{single: []; multiple: []}
  }

  // optional panel of normals file
  pon_file = file(params.mutect2ponFile)

  scatter_count = Channel.from(params.scatter_count)
//   padding = params.readLength + 100

  MiXMHC2PRED   = ( params.MiXMHC2PRED != "" ) ? file(params.MiXMHC2PRED) : ""

  // check HLAHD & OptiType
  have_HLAHD = false
  run_OptiType = (params.disable_OptiType) ? false : true

  if (params.HLAHD_DIR != "") {
      HLAHD = file(params.HLAHD_DIR + "/bin/hlahd.sh")
      if (checkToolAvailable(HLAHD, "exists", "warn")) {
          HLAHD_DIR  = file(params.HLAHD_DIR)
          HLAHD_PATH = HLAHD_DIR + "/bin"
          if(params.HLAHD_module != "") {
              if (checkToolAvailable("bowtie2", "inPath", "warn", module=params.HLAHD_module)) {
                  have_HLAHD = true
              }
          } else {
              if (checkToolAvailable("bowtie2", "inPath", "warn")) {
                  have_HLAHD = true
              }
          }
      }
  }
  if (! have_HLAHD && run_OptiType) {
      log.warn "WARNING: HLAHD not available - can not predict Class II neoepitopes"
  } else if (! have_HLAHD && ! run_OptiType && use_custom_hlas) {
      log.warn "WARNING: HLAHD not available and OptiType disabled - using only user supplied HLA types"
  } else if (! have_HLAHD && ! run_OptiType && ! use_custom_hlas) {
      exit 1, "ERROR: HLAHD not available and OptiType disabled - can not predict HLA types"
  }

  // check if all tools are installed when not running conda or singularity
  have_vep = false
  if (! workflow.profile.contains('conda') && ! workflow.profile.contains('singularity')) {
      def execTools = ["fastqc", "fastp", "bwa", "samtools", "sambamba", "gatk", "vep", "bam-readcount",
                       "perl", "bgzip", "tabix", "bcftools", "yara_mapper", "python", "cnvkit.py",
                       "OptiTypePipeline.py", "alleleCounter", "freec", "Rscript", "java", "multiqc",
                       "sequenza-utils"]

      for (tool in execTools) {
          checkToolAvailable(tool, "inPath", "error")
      }

      VARSCAN = "java -jar " + file(params.VARSCAN)
      have_vep = true
  } else {
      VARSCAN = "varscan "
  }

  // check if we have mutect1 installed
  have_Mutect1 = false
  if (params.MUTECT1 != "" && file(params.MUTECT1) && params.JAVA7 != "" && file(params.JAVA7)) {
      if(checkToolAvailable(params.JAVA7, "exists", "warn") && checkToolAvailable(params.MUTECT1, "exists", "warn")) {
          JAVA7 = file(params.JAVA7)
          MUTECT1 = file(params.MUTECT1)
          have_Mutect1 = true
      }
  }

  // check if we have GATK3 installed
  have_GATK3 = false
  if (params.GATK3 != "" && file(params.GATK3) && params.JAVA8 != "" && file(params.JAVA8) && ! workflow.profile.contains('conda') && ! workflow.profile.contains('singularity')) {
      if(checkToolAvailable(params.JAVA8, "inPath", "warn") && checkToolAvailable(params.GATK3, "exists", "warn")) {
          JAVA8 = file(params.JAVA8)
          GATK3 = file(params.GATK3)
          have_GATK3 = true
      }
  } else if (workflow.profile.contains('singularity')) {
      JAVA8 = "java"
      GATK3 = "/usr/local/opt/gatk-3.8/GenomeAnalysisTK.jar"
      have_GATK3 = true
  } else if (workflow.profile.contains('conda')) {
      JAVA8 = "java"
      GATK3 = "\$CONDA_PREFIX/opt/gatk-3.8/GenomeAnalysisTK.jar"
      have_GATK3 = true
  }


  // check MIXCR licence
  if (params.TCR && params.MIXCR_lic != "") {
      checkToolAvailable(params.MIXCR_lic, "exists", "error")
  } else if (params.TCR && params.MIXCR_lic == "") {
      exit 1, "ERROR: no MiXCR license file specified, please provide a MiXCR license file in params.config or by using the --MIXCR_lic option.\nIf you do not have a MiXCR license you may:\n\ta) run nextNEOpi with --TCR false\n\tb) request one at https://licensing.milaboratories.com"
  }

  /*
    ________________________________________________________________________________

                                P R O C E S S E S
    ________________________________________________________________________________


    TODO: 更改 into 时更换的变量名
      1. MarkDuplicates_out_ch
      2. ScatteredIntervalListToBed_out_ch
      3. RegionsBedToTabix_out_ch
      4. SplitIntervals_out_ch(SplitIntervals_out_xxx_ch)
      5. preprocessIntervalList_out_ch
      6. RegionsBedToIntervalList_out_ch
      7. RegionsBedToTabix_out_ch
      8. BaseRecalGATK4_out_ch/GatherRecalBamFiles_out_IndelRealignerIntervals_ch0
  */
  gatk4_chck_file = file(baseDir + "/assets/.gatk4_install_ok.chck")
  gatk4_chck_file.append()
  nextNEOpiENV_setup_ch0 = Channel.value("OK")


  // meta 处理
  if (bamInput) {
    check_PE(batch_raw_data_ch, nextNEOpiENV_setup_ch0)
    bam_ch = check_PE.out
    check_seqLib_ch = check_PE.out
    seqLibTypes_ok = Channel.value(check_seqLibTypes_ok(check_seqLib_ch))
    bam2fastq(bam_ch.combine(seqLibTypes_ok), nextNEOpiENV_setup_ch0)
    bam_fastq_ch = bam2fastq.out
  } else {
    bam_fastq_ch = channel.empty()
  }


  // TODO: fastq_ch
  fastq_multi_ch = fastq_ch.multiple.map{
        meta, r -> {
            def r1 = []
            def r2 = []
            r.eachWithIndex{ v, ix -> ( ix & 1 ? r2 : r1 ) << v } // TODO: 此处的 << 要进行处理
            tuple( meta, r1.sort(), r2.sort())
        }
  }
  merge_fastq(fastq_multi_ch)
  merged_fastq_ch = merge_fastq.out


  // here we have our final raw fastq files:
  // merged if multi lane runs are used
  // bam2fastq if bam files are used
  raw_reads_ch = merged_fastq_ch.mix(fastq_ch.single, bam_fastq_ch)
  fastqc_reads_ch = merged_fastq_ch.mix(fastq_ch.single, bam_fastq_ch)


  // Common region files preparation for faster processing
  if (params.WES) {
    ref_region_ch = Channel.value([reference.RefDict, reference.RegionsBed])
    regions_bed_to_interval_list(ref_region_ch, nextNEOpiENV_setup_ch0)
    RegionsBedToIntervalList_out_ch = regions_bed_to_interval_list.out

    ref_baits_ch = Channel.value([reference.RefDict, reference.BaitsBed])
    baits_bed_to_interval_list(ref_baits_ch, nextNEOpiENV_setup_ch0)
    BaitsBedToIntervalList_out_ch0 = baits_bed_to_interval_list.out
  } else {
    RegionsBedToIntervalList_out_ch = Channel.empty()
    BaitsBedToIntervalList_out_ch0 = Channel.empty()
  }


  ref_fasta_ch = Channel.value([reference.RefFasta, reference.RefIdx, reference.RefDict])
  preprocess_interval_list(ref_fasta_ch, RegionsBedToIntervalList_out_ch, nextNEOpiENV_setup_ch0)
  preprocessIntervalList_out_ch = preprocess_interval_list.out


  // Splitting interval file in 20(default) files for scattering Mutect2
  split_intervals(ref_fasta_ch, preprocessIntervalList_out_ch, scatter_count, nextNEOpiENV_setup_ch0) // TODO：处理 scatter_count
  SplitIntervals_out_ch = split_intervals.out.path_list
  SplitIntervals_out_ch_Name = split_intervals.out.interval_name


  // convert padded interval list to Bed file (used by varscan)
  // generate a padded tabix indexed region BED file for strelka
  interval_list_to_bed(preprocessIntervalList_out_ch, nextNEOpiENV_setup_ch0)
  RegionsBedToTabix_out_ch = interval_list_to_bed.out


  // convert scattered padded interval list to Bed file (used by varscan)
  scattered_interval_list_to_bed(SplitIntervals_out_ch_Name.combine(SplitIntervals_out_ch.flatten()), nextNEOpiENV_setup_ch0)
  ScatteredIntervalListToBed_out_ch = scattered_interval_list_to_bed.out


  // FastQC
  fast_qc(fastqc_reads_ch)
  ch_fastqc = fast_qc.out // multiQC


  // adapter trimming
  //
  // We first check if we need to trim DNA and RNA or only one of them.
  // If it is only one of them we need to combine the trimmed and raw
  // reads again

  def trim_adapters = false

  if (params.trim_adapters || params.trim_adapters_RNAseq) {
      trim_adapters = true
      reads_to_trim = raw_reads_ch

      if (params.trim_adapters && params.trim_adapters_RNAseq) {
          reads_to_trim_ch = raw_reads_ch
          reads_to_keep_ch = Channel.empty()
      }

      if (params.trim_adapters && ! params.trim_adapters_RNAseq) {
          raw_reads_ch.branch {
              DNA: it[0].sampleType != "tumor_RNA"
              RNA: it[0].sampleType == "tumor_RNA"
          }
          .set{ raw_reads_ch }

          reads_to_trim_ch = raw_reads_ch.DNA
          reads_to_keep_ch = raw_reads_ch.RNA
      }

      if (! params.trim_adapters && params.trim_adapters_RNAseq) {
          raw_reads_ch.branch {
              DNA: it[0].sampleType != "tumor_RNA"
              RNA: it[0].sampleType == "tumor_RNA"
          }
          .set{ raw_reads_ch }

          reads_to_trim_ch = raw_reads_ch.RNA
          reads_to_keep_ch = raw_reads_ch.DNA
      }
  }

  if (trim_adapters) {
    fastp(reads_to_trim_ch)
    reads_trimmed_ch = fastp.out.reads
    fastqc_reads_trimmed_ch = fastp.out.reads
    ch_fastp = fastp.out.json

    // FastQC after adapter trimming
    fast_qc_trimmed(fastqc_reads_trimmed_ch, )
    ch_fastqc_trimmed = fast_qc_trimmed.out

    // combine trimmed reads ch with reads channel of reads that
    // did not need trimming
    reads_ch = reads_trimmed_ch.mix(reads_to_keep_ch)
  } else {
    // no adapter trimming
    ch_fastqc_trimmed = Channel.empty()
    reads_ch = raw_reads_ch
    ch_fastp = Channel.empty()
  }

  // get DNA/RNA reads
  reads_ch.branch {
      DNA: it[0].sampleType != "tumor_RNA"
      RNA: it[0].sampleType == "tumor_RNA"
  }
  .set{ reads_ch }

  reads_BAM_ch = reads_ch.DNA
  reads_uBAM_ch = reads_ch.DNA
  reads_mixcr_DNA_ch = reads_ch.DNA
  dummy_ch = reads_ch.DNA

  reads_tumor_optitype_ch = reads_ch.RNA
  reads_tumor_hlahd_RNA_ch = reads_ch.RNA
  reads_tumor_neofuse_ch = reads_ch.RNA
  reads_tumor_mixcr_RNA_ch = reads_ch.RNA

  reads_mixcr_ch = reads_mixcr_DNA_ch.mix(reads_tumor_mixcr_RNA_ch)

  // setup dummy RNA channels
  no_RNA_ch = dummy_ch.filter{
      it[0].have_RNA == false
  }.map{
      it ->
          return [it[0], []]
  }


  /////// start processing reads ///////

  // make uBAM
  make_ubam(reads_uBAM_ch, nextNEOpiENV_setup_ch0)
  uBAM_out_ch0 = make_ubam.out


  // Aligning reads to reference, sort and index; create BAMs
  ref_bwa_ch = Channel.value([reference.RefFasta, reference.RefIdx, reference.RefDict, reference.BwaRef])
  bwa(reads_BAM_ch, ref_bwa_ch, nextNEOpiENV_setup_ch0)
  BWA_out_ch0 = bwa.out


  // merge alinged BAM and uBAM
  merge_ubam_bam(BWA_out_ch0.join(uBAM_out_ch0, by: [0]), ref_fasta_ch, nextNEOpiENV_setup_ch0)
  uBAM_BAM_out_ch = merge_ubam_bam.out


  // Mark duplicates with sambamba
  mark_duplicates(uBAM_BAM_out_ch, ref_fasta_ch, nextNEOpiENV_setup_ch0)
  MarkDuplicates_out_ch = mark_duplicates.out
  MarkDuplicates_out_ch3 = mark_duplicates.out
  MarkDuplicates_out_ch4 = mark_duplicates.out


  // prepare channel for mhc_extract -> hlad-hd, optitype
  MarkDuplicates_out_ch3.filter {
                                      it[0].sampleType == "tumor_DNA"
                                  }
                                  .set { MarkDuplicatesTumor_out_ch0 }

  // spilt T/N and remove differing/unused info from meta for joining
  // this prepares T/N channel for CNVkit

  MarkDuplicates_out_ch4.branch {
          meta_ori, bam ->
              def meta = meta_ori.clone()
              tumor : meta.sampleType == "tumor_DNA"
                  meta.remove('sampleType')
                  return [meta, bam]

              normal: meta.sampleType == "normal_DNA"
                  meta.remove('sampleType')
                  return [meta, bam]
  }.set{ MarkDuplicates_out_ch4 }

  MarkDuplicates_out_CNVkit_ch0 = MarkDuplicates_out_ch4.tumor.join(MarkDuplicates_out_ch4.normal, by:[0])


  if (params.WES) {
    // Generate HS metrics using picard
    ref_alig_ch = Channel.value([reference.RefFasta, reference.RefIdx])
    alignment_metrics(
        MarkDuplicates_out_ch, ref_alig_ch, BaitsBedToIntervalList_out_ch0, 
        RegionsBedToIntervalList_out_ch, nextNEOpiENV_setup_ch0
    )
    alignmentMetrics_ch = alignment_metrics.out
  } else {
    // bogus channel for multiqc
    alignmentMetrics_ch = Channel.empty()
  }


  /*
   BaseRecalibrator (GATK4): generates recalibration table for Base Quality Score
   Recalibration (BQSR)
   ApplyBQSR (GATK4): apply BQSR table to reads
  */
  bqsr_db_ch = Channel.value(
    [ database.MillsGold,
      database.MillsGoldIdx,
      database.DBSNP,
      database.DBSNPIdx,
      database.KnownIndels,
      database.KnownIndelsIdx ]
  )
  scatter_base_recal_gatk4(
    MarkDuplicates_out_ch.combine(SplitIntervals_out_ch.flatten()),
    ref_fasta_ch, bqsr_db_ch, nextNEOpiENV_setup_ch0
  )
  scatterBaseRecalGATK4_out_ch0 = scatter_base_recal_gatk4.out


  // gather scattered bqsr tables
  gather_gatk4_scsattered_bqsr_tables(
    scatterBaseRecalGATK4_out_ch0.groupTuple(by: [0]), nextNEOpiENV_setup_ch0
  )
  gatherBQSRtables_out_ch0 = gather_gatk4_scsattered_bqsr_tables.out


  // ApplyBQSR (GATK4): apply BQSR table to reads
  scatter_gatk4_apply_bqsrs(
    MarkDuplicates_out_ch.join(gatherBQSRtables_out_ch0, by: [0])
    .combine(
        SplitIntervals_out_ch.flatten()
    ),
    ref_fasta_ch, bqsr_db_ch, nextNEOpiENV_setup_ch0
  )
  scatterGATK4applyBQSRS_out_GatherRecalBamFiles_ch0 = scatter_gatk4_apply_bqsrs.out


  gather_recal_bam_files(
    scatterGATK4applyBQSRS_out_GatherRecalBamFiles_ch0
    .toSortedList({a, b -> a[1][0].baseName <=> b[1][0].baseName})
    .flatten().collate(3).groupTuple(by: [0]),
    nextNEOpiENV_setup_ch0
  )
  BaseRecalGATK4_out_ch = gather_recal_bam_files.out


  // GetPileupSummaries (GATK4): tabulates pileup metrics for inferring contamination
  // get_pileup
  // get_pileup

}