# Insertion-Deletion Calling Workflow for Roddy

A Platypus-based insertion/deletion-detection workflow with extensive quality control additions for the workflow management system [Roddy](https://github.com/eilslabs/Roddy).

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

# Run flags / switches

Switch                      Default Description
runIndelAnnotation          true    Run the annotation step or stop the workflow before it.
runIndelDeepAnnotation      true    Run the deep annotation step or stop the workflow before it.
runIndelVCFFilter           true    Run the filter step or stop the workflow before it.
runTinda                    true    Check for sample swaps with TiNDA.

# Changelist

* Version update to 1.2.177-602

- Minor: Adding the local control annotation (WGS + WES) via the deepAnnotation pipe. These annotations will be used downstream for filtering variants from the reduced control coverage settings in the WGS samples.

* Version update to 1.2.177-601

- Major: Fix the control-tumor column swap bug

* Version update to 1.2.177

- Roddy 3.0 support
- 1.2.177 is equivalent to 1.0.176-9
- 1.2.177-1 is equivalent to 1.0.176-10.

* Version update to 1.0.176-9

- Fix of a bug affecting versions 1.0.176 to 1.0.176-8 that results in higher false positive rate.

* Version update to 1.0.176

- SNVs calling made default. 
- Swapchecker - checks for tumor/control swap from the same PID. 
- TiNDA - Tumor in normal detection analysis, using Canopy's EM-clustering algorithm

* Version update to 1.0.168

- Further checks for the platypus indel calling step are introduced. A zgrep will be performed, together with a linecount to
  see if there are any faulty lines in the raw vcf file.

* Version update to 1.0.161

* Version update to 1.0.157

- Move the indel calling workflow to its own plugin

* Version update to 1.0.131

- Change workflow class to override another execute method. This makes the workflow a bit cleaner.