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

* Local controls:
  * VCF files containing the frequency of the variants collected from the control samples sequenced and aligned locally using the same/similar workflows.
  * Separate files for WES and WGS control samples.
  * INFO column should contain the AF field.

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

* Version update to 3.1.0 (2.6.0-deprecated) 

  - Bugfix: Fix of column swap bug introduced before 1.0.167. Output VCF with swapped control and tumor genotype columns if they are in the 11th and 10th column respectively.
  - Crash the workflow if genotype column names could not be verified through BAM SM tags
  - Added `--skip_order_tag` to skip the above if users manually verified the genotype column order
  - Crash the workflow if more than two samples are present in the raw VCF file. This could be due to multiple RG tags.

* Version update to 3.0.0 (2.5.0-deprecated)

  - Added a local control generated from ~1k WES samples
  
* Version update to 2.4.3

  - Bugfix: `platypusIndelAnnotation.sh` now checks output of annovar execution for error code != 0. Without this the workflow will continue even if annovar throws an exception.

* Version update to 2.4.2

  - Fixed missed annotations due to unsorted VCF

* Version update to 2.4.1-1 (ReleaseBranch_2.4.1)

  - Bugfix: Fix of column swap bug introduced before 1.0.167. Output VCF with swapped control and tumor genotype columns if they are in the 11th and 10th column respectively.

* Version update to 2.4.1

  - Create __non-empty__ sample swap JSON even if less than 50 germline variants
  - More error-robust IO in sample swap TINDA Perl script

* Version update to 2.4.0

  - Bugfix: Casting error in canopy-based clustering R code
  - Bugfix: Create warning PDF with text if number of to-be-merged PDFs is too large or zero, plus fix one-off bug
  - Exit 0 instead of != 0, if less then 50 germline variants
 
* Version update 2.3.0

  - Added `isNoControlWorkflow` variable and make FILTER_database work with it
  - Removed usage of ExAC from filtering, gnomAD includes ExAC
  - Report only exact matches for database annotation

* Version update to 2.2.2 (ReleaseBranch_2.2)

  - Bugfix: Fix of column swap bug introduced before 1.0.167. Output VCF with swapped control and tumor genotype columns if they are in the 11th and 10th column respectively.

* Version update to 2.2.1 
 
  - Stability improvement Perl to protect against I/O errors
  - Write text PDf upon empty and too many indels

* Version update to 2.2.0

  - Upgrade from [COWorkflowsBasePlugin](https://github.com/DKFZ-ODCF/COWorkflowsBasePlugin) 1.1.0 to 1.4.1

* Version update to 2.1.0-2

  - Fixed FILENAME_ parameter
  - Executability check for `REFERENCE_GENOME` variable (file accessible from submission host)
  - Fixed reported version.
  - Code cleanup in `annotate_vcf.pl`

* Version update to 2.1.0-1

  - Added gnomAD exomes
  - Added local controls
  - Adding python script for parsing MNPs
  - Check REF_genome for BAM files
  - Check BAM file readability

* Version update to 2.0.0-2

  - Restricting screenshot generation to 100
  
* Version update to 2.0.0-1

  - README update

* Version update to 2.0.0-101 (ReleaseBranch_2.0.0-1)

  - Bugfix: Fix of column swap bug introduced before 1.0.167. Output VCF with swapped control and tumor genotype columns if they are in the 11th and 10th column respectively.
  
* Version update to 2.0.0

  - TiNDA workflow was updated
    - Two local controls are removed, and a new local control created from ~3000 WGS platypus variant calls were added.
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

* Version update to 1.2.177-8

* Version update to 1.2.177-7

* Version update to 1.2.177-601 (ReleaseBranch_1.2.177-6)

  - Bugfix: Fix of column swap bug introduced before 1.0.167. Output VCF with swapped control and tumor genotype columns if they are in the 11th and 10th column respectively.

* Version update to 1.2.177-6

  - removed a `umask` from `checkSampleSwap_TiN.sh`

* Version update to 1.2.177-5

  - Platypus update 0.8.1 to 0.8.1.1 (this only fixes the reported Platypus version number in the VCF)
  - Added conda environment
  - Reordered cvalues in plugin XML configuration
  - Neutral refactoring in `indelCalling.sh`
  - Bugfix: JSON syntax for "file"

* Version update to 1.2.177-4 (includes development versions 1.2.177-1 to -3)

  - Roddy 3.0 support
  - Updated COWorkflows base plugin from 1.2.59 to 1.2.76
  - Fixed reported version
  - Bugfix in `indel_extractor_v1.pl` concerning ncRNA_exonic and ncRNA_splicing annotations for somatic calls.
  - Added `runTindo` cvalue

* Version update to 1.0.176-10 (=1.2.177-1)

  - Added right and bottom border features to `TiN_canopyBasedClustering.R`
  - Version is not backwards compatible

* Version update to 1.0.176-9 (=1.0.177)

  - Fix of a bug affecting versions 1.0.176 to 1.0.176-8 that results in higher false positive rate.

* Version update to 1.0.176

  - SNVs calling made default.
  - Swapchecker - checks for tumor/control swap from the same PID.
  - TiNDA - Tumor in normal detection analysis, using Canopy's EM-clustering algorithm.

* Version update to 1.0.168

  - Further checks for the platypus indel calling step are introduced. A zgrep will be performed, together with a linecount to see if there are any faulty lines in the raw vcf file.
  
* Version update to 1.0.167

  - Note that commit e8d8fc27 _introduced_ a bug. Output VCF with swapped control and tumor genotype columns if they are in the 11th and 10th column respectively.

* Version update to 1.0.161

* Version update to 1.0.157

  - Move the indel calling workflow to its own plugin.

* Version update to 1.0.131

  - Change workflow class to override another execute method. This makes the workflow a bit cleaner.
