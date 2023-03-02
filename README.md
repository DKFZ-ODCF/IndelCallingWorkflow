# Insertion-Deletion Calling Workflow for Roddy

A Platypus-based insertion/deletion-detection workflow with extensive quality control additions for the workflow management system [Roddy](https://github.com/eilslabs/Roddy). The workflow is only suited for human data (hg37/hg19; some later versions also hg38), because of the big role annotations play in this workflow.

> <table><tr><td><a href="https://www.denbi.de/"><img src="docs/images/denbi.png" alt="de.NBI logo" width="300" align="left"></a></td><td><strong>Your opinion matters!</strong> The development of this workflow is supported by the <a href="https://www.denbi.de/">German Network for Bioinformatic Infrastructure (de.NBI)</a>. By completing <a href="https://www.surveymonkey.de/r/denbi-service?sc=hd-hub&tool=IndelCallingWorkflow">this very short (30-60 seconds) survey</a> you support our efforts to improve this tool.</td></tr></table>

## Citing

The Indel workflow was in the pan-cancer analysis of whole genomes (PCAWG) and can be cited with the following publication:

* Pan-cancer analysis of whole genomes.<br>
  The ICGC/TCGA Pan-Cancer Analysis of Whole Genomes Consortium.<br>
  Nature volume 578, pages 82–93 (2020).<br>
  DOI [10.1038/s41586-020-1969-6](https://doi.org/10.1038/s41586-020-1969-6)

Containers are available in [Dockstore](https://dockstore.org/containers/quay.io/pancancer/pcawg-dkfz-workflow).

Furthermore, the workflow is regularly used at the DKFZ throug the automation system [One-Touch Pipeline](https://gitlab.com/one-touch-pipeline/otp):

* OTP: An automatized system for managing and processing NGS data.<br>
  Eva Reisinger, Lena Genthner, Jules Kerssemakers, Philip Kensche, Stefan Borufka, Alke Jugold, Andreas Kling, Manuel Prinz, Ingrid Scholz, Gideon Zipprich, Roland Eils, Christian Lawerenz, Jürgen Eils.<br>
   Journal of Biotechnology, volume 261, pages 53-62 (2017).<br>
   DOI: [10.1016/j.jbiotec.2017.08.006](https://doi.org/10.1016/j.jbiotec.2017.08.006)

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

## General Remarks 

Changes in the TiNDA and no-control functionalities do not trigger a major version bump.

The currently supported branches are 3 (master), 2.4.1, 2.2, 2.0.0-1, and 1.2.177. Unfortunately, the results from these branches are not compatible with each other. For larger studies you may want to use only versions that produce comparable results.

Note that only on `master` new features are implemented, so the other branches are only for bugfixes.

## Change-Log

* Version update to 3.1.1

  - Patch (Bugfix): The nocontrol workflow is exempted from the tumor & control column swap introduced in 3.1.0.

* Version update to 3.1.0 (2.6.0-deprecated) 

* Version update to 3.1.0 (2.6.0-deprecated, master)

  - Minor: Bugfix: Fix of column swap [bug](#column-bug) introduced before 1.0.167. Output VCF with swapped control and tumor genotype columns if they are in the 11th and 10th column respectively.
  - Minor: Added `--skip_order_tag` to skip the above, if users manually verified the genotype column order
  - Patch: Crash the workflow if genotype column names could not be verified through BAM SM tags
  - Patch: Crash the workflow if more than two samples are present in the raw VCF file. This could be due to multiple RG tags.
  - Deprecated version-tag 2.6.0 because of the major-level change in 2.5.0-deprecated/3.0.0

* Version update to 2.4.1-1 (ReleaseBranch_2.4.1)

  - Bugfix: Fix of column swap [bug](#column-bug) introduced before 1.0.167. Output VCF with swapped control and tumor genotype columns if they are in the 11th and 10th column respectively.

* Version update to 2.2.2 (ReleaseBranch_2.2)

  - Bugfix: Fix of column swap [bug](#column-bug) introduced before 1.0.167. Output VCF with swapped control and tumor genotype columns if they are in the 11th and 10th column respectively.

* Version update to 2.0.0-101 (ReleaseBranch_2.0.0-1)

  - Bugfix: Fix of column swap [bug](#column-bug) introduced before 1.0.167. Output VCF with swapped control and tumor genotype columns if they are in the 11th and 10th column respectively.

* Version update to 1.2.177-601 (ReleaseBranch_1.2.177-6)

  - Bugfix: Fix of column swap [bug](#column-bug) introduced before 1.0.167. Output VCF with swapped control and tumor genotype columns if they are in the 11th and 10th column respectively.


The following are older versions of the workflow on which no further development will take place.


* Version update to 3.0.0 (2.5.0-deprecated) ([Column bug](#column-bug))

  - Major: Added a local control generated from ~1k WES samples
  - Deprecated version-tag 2.5.0 because this is actually a major-level change.

* Version update to 2.4.3 ([Column bug](#column-bug))

  - Bugfix: `platypusIndelAnnotation.sh` now checks output of annovar execution for error code != 0. Without this the workflow will continue even if annovar throws an exception.

* Version update to 2.4.2 ([Column bug](#column-bug))

  - Fixed missed annotations due to unsorted VCF (sort multi-SNP after variant normalization; only TiNDA affected, but not significantly)
  
* Version update to 2.4.1 ([Column bug](#column-bug))

  - Create __non-empty__ sample swap JSON even if less than 50 germline variants
  - More error-robust IO in sample swap TINDA Perl script

* Version update to 2.4.0 ([Column bug](#column-bug))

  - Bugfix: Casting error in canopy-based clustering R code
  - Bugfix: Create warning PDF with text if number of to-be-merged PDFs is too large or zero, plus fix one-off bug
  - Exit 0 instead of != 0, if less then 50 germline variants

* Version update 2.3.0 ([Column bug](#column-bug))

  - Added `isNoControlWorkflow` variable and make FILTER_database work with it
  - Removed usage of ExAC from filtering, gnomAD includes ExAC (only no-control affected)
  - Report only exact matches for database annotation (only no-control affected)
  
* Version update to 2.2.1 ([Column bug](#column-bug))

  - Stability improvement Perl to protect against I/O errors
  - Write text PDF upon empty and too many indels

* Version update to 2.2.0 ([Column bug](#column-bug))

  - Upgrade from [COWorkflowsBasePlugin](https://github.com/DKFZ-ODCF/COWorkflowsBasePlugin) 1.1.0 to 1.4.1

* Version update to 2.1.0-2 ([Column bug](#column-bug))

  - Fixed FILENAME_ parameter
  - Executability check for `REFERENCE_GENOME` variable (file accessible from submission host)
  - Fixed reported version.
  - Code cleanup in `annotate_vcf.pl`

* Version update to 2.1.0-1 ([Column bug](#column-bug))

  - Added gnomAD exomes
  - Added local controls
  - Adding python script for parsing MNPs
  - Check REF_genome for BAM files
  - Check BAM file readability

* Version update to 2.0.0-2 ([Column bug](#column-bug))

  - Restricting screenshot generation to 100

* Version update to 2.0.0-1 ([Column bug](#column-bug))

  - README update
  
* Version update to 2.0.0 ([Column bug](#column-bug))

  - TiNDA workflow was updated
    - Two local controls are removed, and a new local control created from ~3000 WGS platypus variant calls were added.
    - gnomAD v2.1 exomes and genomes are added
    - Variants with MAF above 0.01 in any one of 3 of the background data sets (gnomAD exomes, gnomAD genomes or local control) were considered as common SNPs are were removed.
    - Rare germline variants and somatic variants are annotated with ANNOVAR using gencode v19 model
    - Temporary TiNDA files are cleaned up
    - TiNDA will stop if there are less than 50 rare germline variants

* Version update to 1.3.0 ([Column bug](#column-bug))

  - Added `tumorSample` and `controlSample` variables to job environments. These can also be used in output file paths.

* Version update to 1.2.178 ([Column bug](#column-bug))

  - Including gnomAD exomes and genomes for nocontrol workflow filtering.
  - Adding gnomAD files for no-control workflow.
  - Updated from COWorkflowsBasePlugin:1.0.3 to COWorkflowsBasePlugin:1.1.0.

* Version update to 1.2.177-8 (=1.2.182) ([Column bug](#column-bug))

  - major: Changed `GERMLINE_AVAILABLE` with allowed values 0 and 1 to `isNoControlWorkflow` with allowed values "false" and "true".
  - minor: Update to from COWorkflows:1.2.59 to COWorkflowsBasePlugin:1.0.3.
  - minor: Added writing TiN classification to the VCF file.
  - Note that the version may report itself as 1.2.182 according to the `versionInfo.txt`.

* Version update to 1.2.177-7 ([Column bug](#column-bug))

  - patch: Fix problem with dealing with empty result files in `indelCalling.sh`.
  
* Version update to 1.2.177-6 ([Column bug](#column-bug))

  - Removed a `umask` from `checkSampleSwap_TiN.sh`.

* Version update to 1.2.177-5 ([Column bug](#column-bug))

  - Platypus update 0.8.1 to 0.8.1.1 (this only fixes the reported Platypus version number in the VCF)
  - Added conda environment
  - Reordered cvalues in plugin XML configuration
  - Neutral refactoring in `indelCalling.sh`
  - Bugfix: JSON syntax for "file"

* Version update to 1.2.177-4 (includes development versions 1.2.177-1 to -3) ([Column bug](#column-bug))

  - Roddy 3.0 support
  - Updated COWorkflows base plugin from 1.2.59 to 1.2.76
  - Fixed reported version
  - Bugfix in `indel_extractor_v1.pl` concerning ncRNA_exonic and ncRNA_splicing annotations for somatic calls.
  - Added `runTinda` cvalue

* Version update to 1.0.177-2 ([Column bug](#column-bug))

  - Bugfix ncRNA file
  - Allow toggling on/off of Tinda

* Version update to 1.0.176-10 (=1.0.177-1) ([Column bug](#column-bug))

  - Major: Added a method based on maximum Control_AF to remove unusual cluster.
  - Major: Increased the number of iterations to 100 (from 50) in the EM step. TiNDA rescue clusters will not be same.
  - Minor: Added a polygon shape to the TiNDA plot. Plots will look different.
  - Minor: Parametrised the right and left threshold that defines the border for TiNDA cluster selection, but by default use the same values as before.
  - Patch: Fixed the bugs in the bash sanity checks.

* Version update to 1.0.176-9 (=1.0.177) ([Column bug](#column-bug))

  - Fix of a bug affecting versions 1.0.176 to 1.0.176-8 that results in higher false positive rate.

> NOTE: Refrain from using versions 1.0.176 through 1.0.176-8. They are affected by a bug that increases the false positive rate.

* Version update to 1.0.176 ([Column bug](#column-bug))

  - SNVs calling made default.
  - Swapchecker - checks for tumor/control swap from the same PID.
  - TiNDA - Tumor in normal detection analysis, using Canopy's EM-clustering algorithm.

* Version update to 1.0.168 ([Column bug](#column-bug))

  - Further checks for the platypus indel calling step are introduced. A zgrep will be performed, together with a linecount to see if there are any faulty lines in the raw vcf file.

* Version update to 1.0.167 (<a id="column-bug"></a>Column bug)

  - Background: Platypus always orders the columns alpha-numerically. In our data, usually sample ID for the control is alpha-numerically smaller than that for the tumor, e.g., control/normal/blood < tumor/metastasis. Only for the sample pairs where tumor ID comes before the control (e.g. 'cell_line' < 'control'), the tumor genotype columns come before the control column in the raw VCF. This column swap was accounted for in the workflow code while detecting the somatic variants: if needed, the columns were swapped. 
    Bug: The commit e8d8fc27, which replaced the Perl script with the Python script, introduced a bug that replaced the names of the columns in the 10th & 11th columns with the unsorted labels, but left the data in the order created by Platypus.
    > NOTE: This bug has no influence on the list of reported indels, and all somatic indels are correct. Only the genotype data are in the wrong columns.

* Version update to 1.0.161

* Version update to 1.0.157

  - Move the indel calling workflow to its own plugin.

* Version update to 1.0.131

  - Change workflow class to override another execute method. This makes the workflow a bit cleaner.
