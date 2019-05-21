# Insertion-Deletion Calling Workflow for Roddy

A Platypus-based insertion/deletion-detection workflow with extensive quality control additions for the workflow management system [Roddy](https://github.com/eilslabs/Roddy).

> <table><tr><td><a href="https://www.denbi.de/"><img src="docs/images/denbi.png" alt="de.NBI logo" width="300" align="left"></a></td><td><strong>Your opinion matters!</strong> The development of this workflow is supported by the <a href="https://www.denbi.de/">German Network for Bioinformatic Infrastructure (de.NBI)</a>. By completing <a href="https://www.surveymonkey.de/r/denbi-service?sc=hd-hub&tool=IndelCallingWorkflow">this very short (30-60 seconds) survey</a> you support our efforts to improve this tool.</td></tr></table>

## Software Requirements

Please refer to the Roddy [website](https://github.com/eilslabs/Roddy) for instructions on how to install Roddy.

Most bioinformatic software required for this workflow is available in Conda, however, you additionally need Pypy -- preferentially version >= 5. You can install it via your operating system's package manager. Alternatively, you can install a portable Pypy from [squeaky-pl](https://github.com/squeaky-pl/portable-pypy), in particular if the OS version is too old (currently we use pypy 5.8).

### Conda

The workflow contains a description of a [Conda](https://conda.io/docs/) environment. A number of Conda packages from [BioConda](https://bioconda.github.io/index.html) are required. You should set up the Conda environment at a centralized position available from all compute hosts. 

First install the BioConda channels:
```
conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```

Then install the environment

```
conda env create -n IndelCallingWorkflow -f $PATH_TO_PLUGIN_DIRECTORY/resources/analysisTools/indelCallingWorkflow/environments/conda.yml
```

The name of the Conda environment is arbitrary but needs to be consistent with the `condaEnvironmentName` variable. The default for that variable is set in `resources/configurationFiles/analysisIndelCalling.xml`.

## Data Requirements

There are quite extensive requirements in annotation etc. data required for the workflow. Please have a look at the file `resources/configurationFiles/analysisIndelCalling.xml`. Note that input all VCF and BED files need to be indexed with tabix.

# Configuration Values

|Switch                    |  Default     | Description
|--------------------------|--------------|-----------------------------------------------|
| runIndelAnnotation       |  true        | Run the annotation step or stop the workflow before it. |
| runIndelDeepAnnotation   |  true        | Run the deep annotation step or stop the workflow before it. |
| runIndelVCFFilter        |  true        | Run the filter step or stop the workflow before it. |
| runTinda                 |  true        | Check for sample swaps with TiNDA. |
| bamfile_list             | empty        | Semicolon-separated list of BAM files, starting with the control's BAM. Each BAM file needs an index file with the same name as the BAM, but ".bai" suffixed. |
| sample_list              | empty        | Semicolon-separated list of sample names in the same order as `bamfile_list` |
| possibleTumorSampleNamePrefixes | "( tumor )" | Bash-array of tumor sample name prefixes. |
| possibleControlSampleNamePrefixes | "( control )" | Bash-array of control sample name prefixes. |
| REFERENCE_GENOME | empty | |
| CHR_SUFFIX | "" | Suffix added to the chromosome names |
| CHR_PREFIX | "" | Prefix added to the chromosome names |
| extractSamplesFromOutputFiles | true | |
| CHROMOSOME_INDICES | empty | Bash-array of chromosome names to which the analysis should be restricted |

Since version 2.2.0 the workflow uses the [COWorkflowsBasePlugin](https://github.com/DKFZ-ODCF/COWorkflowsBasePlugin) 1.4.1+ with an alternative algorithm for extracting sample names from BAM files.

# Example call

TBD

# Changelist

* Version update to 2.2.0

  * Upgrade from [COWorkflowsBasePlugin](https://github.com/DKFZ-ODCF/COWorkflowsBasePlugin) 1.1.0 to 1.4.1
  * Executability check for `REFERENCE_GENOME` variable (file accessible from submission host)

* Version update to 2.1.0-1

  * Added gnomAD exomes
  * Check REF_genome for BAM files
  * Check BAM file readability
  * Added local controls
  * Adding python script for parsing MNPs

* Version update to 2.0.0-2

  * Restricting screenshot generation to 100

* Version update to 2.0.0-1

  * README update

* Version update to 2.0.0

  - TiNDA workflow was updated
   - Two local controls are removed and a new local control created from ~3000 WGS platypus variant calls were added.
   - gnomAD v2.1 exomes and genomes are added
   - Variants with MAF above 0.01 in any one of 3 of the background data sets (gnomAD exomes, gnomAD genomes or local control) were considered as common SNPs are were removed.
   - Rare germline variants and somatic variants are annotated with ANNOVAR using gencode v19 model
   - Temporary TiNDA files are cleaned up
   - TiNDA will stop if there are less than 50 rare germline variants

* Version update to 1.3.0

  - Added `tumorSample` and `controlSample` variables to job environments. These can also be used in output file paths.

* Version update to 1.2.178

  - Including gnomAD exomes and genomes for nocontrol workflow filtering.
  - Adding gnomAD files for no-control workflow.
  - Updating the COWorkflowsBasePlugin.
  - Including gnomAD exomes and genomes for nocontrol workflow filtering

* Version update to 1.2.177

  - Roddy 3.0 support.
  - 1.2.177 is equivalent to 1.0.176-9.
  - 1.2.177-1 is equivalent to 1.0.176-10.

* Version update to 1.0.176-9

  - Fix of a bug affecting versions 1.0.176 to 1.0.176-8 that results in higher false positive rate.

* Version update to 1.0.176

  - SNVs calling made default.
  - Swapchecker - checks for tumor/control swap from the same PID.
  - TiNDA - Tumor in normal detection analysis, using Canopy's EM-clustering algorithm.

* Version update to 1.0.168

  - Further checks for the platypus indel calling step are introduced. A zgrep will be performed, together with a linecount to see if there are any faulty lines in the raw vcf file.

* Version update to 1.0.161

* Version update to 1.0.157

  - Move the indel calling workflow to its own plugin.

* Version update to 1.0.131

  - Change workflow class to override another execute method. This makes the workflow a bit cleaner.
