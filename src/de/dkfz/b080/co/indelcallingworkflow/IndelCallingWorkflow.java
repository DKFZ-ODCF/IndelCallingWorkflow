package de.dkfz.b080.co.indelcallingworkflow;

import de.dkfz.b080.co.common.WorkflowUsingMergedBams;
import de.dkfz.b080.co.files.*;
import de.dkfz.roddy.config.RecursiveOverridableMapContainerForConfigurationValues;
import de.dkfz.roddy.core.ExecutionContext;
import de.dkfz.roddy.knowledge.files.FileObject;
import de.dkfz.roddy.knowledge.files.Tuple2;
import de.dkfz.roddy.tools.LoggerWrapper;

/**
 * Indel calling based on the platypus pipeline.
 */
public class IndelCallingWorkflow extends WorkflowUsingMergedBams {

    private static LoggerWrapper logger = LoggerWrapper.getLogger(IndelCallingWorkflow.class.getName());

    private static VCFFileForIndels extractVCF(FileObject runResult) {
        return (VCFFileForIndels) ((Tuple2)runResult).value0;
    }
    
    @Override
    public boolean execute(ExecutionContext context, BasicBamFile _bamControlMerged, BasicBamFile _bamTumorMerged) {

        boolean runTinda = getFlag("runTinda", true);
        boolean runAnnotation = getFlag("runIndelAnnotation", true);
        boolean runDeepAnnotation = getFlag("runIndelDeepAnnotation", true);
        boolean runFilter = getFlag("runIndelVCFFilter", true);

        RecursiveOverridableMapContainerForConfigurationValues configurationValues = context.getConfiguration().getConfigurationValues();

        TumorBamFile bamTumorMerged = new TumorBamFile(_bamTumorMerged);

        if (isControlWorkflow())
            return executeWithControl(new ControlBamFile(_bamControlMerged), bamTumorMerged, runTinda, runAnnotation, runDeepAnnotation, runFilter);
        else
            return executeWithoutControl(bamTumorMerged, runTinda, runAnnotation, runDeepAnnotation, runFilter);
    }

    private boolean executeWithControl(ControlBamFile bamControlMerged, TumorBamFile bamTumorMerged, boolean runTinda, boolean runAnnotation, boolean runDeepAnnotation, boolean runFilter) {
        VCFFileForIndels rawVCF = (VCFFileForIndels) run("indelCalling", bamTumorMerged, bamControlMerged);

        if (runTinda)
            run("checkSampleSwap", rawVCF, bamTumorMerged, bamControlMerged);

        if (!runAnnotation) return true;
        VCFFileForIndels annotatedVCFFile = extractVCF(run("indelAnnotation", bamTumorMerged, bamControlMerged, rawVCF));

        if (!runDeepAnnotation) return true;
        VCFFileForIndels deepAnnotatedVCFFile = extractVCF(run("indelDeepAnnotation", annotatedVCFFile, "PIPENAME=INDEL_DEEPANNOTATION"));

        if (!runFilter) return true;
        run("indelVcfFilter", bamTumorMerged, bamControlMerged, deepAnnotatedVCFFile);
        return true;
    }

    private boolean executeWithoutControl(TumorBamFile bamTumorMerged, boolean runTinda, boolean runAnnotation, boolean runDeepAnnotation, boolean runFilter) {
        VCFFileForIndels rawVCF = (VCFFileForIndels) run("indelCallingWithoutControl", bamTumorMerged);

        if (runTinda)
            logger.always("Not running Tinda, since no germline.");

        if (!runAnnotation) return true;
        VCFFileForIndels annotatedVCFFile = extractVCF(run("indelAnnotationWithoutControl", bamTumorMerged, rawVCF));

        if (!runDeepAnnotation) return true;
        VCFFileForIndels deepAnnotatedVCFFile = extractVCF(run("indelDeepAnnotation", annotatedVCFFile, "PIPENAME=INDEL_DEEPANNOTATION"));

        if (!runFilter) return true;
        run("indelVcfFilterWithoutControl", bamTumorMerged, deepAnnotatedVCFFile);
        return true;
    }
}
