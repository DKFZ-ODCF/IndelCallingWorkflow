package de.dkfz.b080.co.indelcallingworkflow;

import de.dkfz.b080.co.common.WorkflowUsingMergedBams;
import de.dkfz.b080.co.files.*;
import de.dkfz.roddy.config.RecursiveOverridableMapContainerForConfigurationValues;
import de.dkfz.roddy.core.ExecutionContext;
import de.dkfz.roddy.knowledge.files.Tuple2;

/**
 * Indel calling based on the platypus pipeline.
 */
public class IndelCallingWorkflow extends WorkflowUsingMergedBams {

    @Override
    public boolean execute(ExecutionContext context, BamFile bamControlMerged, BamFile bamTumorMerged) {

        RecursiveOverridableMapContainerForConfigurationValues configurationValues = context.getConfiguration().getConfigurationValues();
        boolean runFilter = configurationValues.getBoolean("runIndelVCFFilter", true);
        boolean runDeepAnnotation = configurationValues.getBoolean("runIndelDeepAnnotation", true);
        boolean runAnnotation = configurationValues.getBoolean("runIndelAnnotation", true);

        VCFFileForIndels rawVCF = bamTumorMerged.callIndels(bamControlMerged);

        if (!runAnnotation) return true;
        VCFFileForIndels vcfFileForIndels = rawVCF; //Use the raw vcf for further processing.
        VCFFileForIndels annotatedVCFFile = vcfFileForIndels.annotate().value0;

        if (!runDeepAnnotation) return true;
        VCFFileForIndels deepAnnotatedVCFFile = annotatedVCFFile.deepAnnotate().value0;

        if (!runFilter) return true;
        deepAnnotatedVCFFile.filter();

        return true;
    }
}
