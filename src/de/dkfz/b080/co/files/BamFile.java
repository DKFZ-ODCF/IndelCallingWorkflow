package de.dkfz.b080.co.files;

import de.dkfz.roddy.knowledge.files.BaseFile;

/**
 * Created by heinold on 15.01.16.
 */
public class BamFile extends BasicBamFile {
    public BamFile(ConstructionHelperForBaseFiles helper) {
        super(helper);
    }

    public BamFile(BaseFile parent) {
        super(parent);
    }

    public BamFile(BasicBamFile parent) {
        super(parent);
    }
}
