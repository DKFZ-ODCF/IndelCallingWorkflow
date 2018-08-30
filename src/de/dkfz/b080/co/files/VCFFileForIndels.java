/*
 * Copyright (c) 2018 German Cancer Research Center (DKFZ).
 *
 * Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/IndelCallingWorkflow).
 *
 */
package de.dkfz.b080.co.files;

import de.dkfz.roddy.knowledge.files.BaseFile;
import de.dkfz.roddy.knowledge.files.Tuple2;
import de.dkfz.roddy.knowledge.methods.GenericMethod;

/**
 *
 * @author michael
 */
public class VCFFileForIndels extends BaseFile {

    public VCFFileForIndels(ConstructionHelperForBaseFiles helper) {
        super(helper);
    }

    public VCFFileForIndels(TumorBamFile parentFile) {
        super(parentFile);
    }

    /**
     * For i.e. somatic indel file
     * @param parentFile
     */
    public VCFFileForIndels(VCFFileForIndels parentFile) {
        super(parentFile);
    }

}
