package de.dkfz.b080.co.files;

import de.dkfz.roddy.core.ExecutionContext;
import de.dkfz.roddy.knowledge.files.BaseFile;
import de.dkfz.roddy.knowledge.files.FileStageSettings;
import de.dkfz.roddy.execution.jobs.JobResult;
import java.io.File;
import java.util.List;

/**
 *
 * @author michael
 */
public class IndelDeepAnnotationFile extends BaseFile {

    public IndelDeepAnnotationFile(File path, ExecutionContext executionContext, JobResult creatingJobsResult, List<BaseFile> parentFiles, FileStageSettings settings) {
        super(path, executionContext, creatingJobsResult, parentFiles, settings);
    }

}
