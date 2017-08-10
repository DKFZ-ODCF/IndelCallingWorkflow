== Description

Some description

== Run flags / switches

Switch                      Default Description
runIndelAnnotation          true    Run the annotation step or stop the workflow before it.
runIndelDeepAnnotation      true    Run the deep annotation step or stop the workflow before it.
runIndelVCFFilter           true    Run the filter step or stop the workflow before it.

== Changelist

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
