/*
 * Copyright (c) 2018 German Cancer Research Center (DKFZ).
 *
 * Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/IndelCallingWorkflow).
 *
 */
package de.dkfz.b080.co.files;

import de.dkfz.roddy.knowledge.files.BaseFile;

/**
 * Created by heinold on 15.01.16.
 */
public class TumorBamFile extends BasicBamFile {
    public TumorBamFile(ConstructionHelperForBaseFiles helper) {
        super(helper);
    }

    public TumorBamFile(BaseFile parent) {
        super(parent);
    }

    public TumorBamFile(BasicBamFile parent) {
        super(parent);
    }
}
