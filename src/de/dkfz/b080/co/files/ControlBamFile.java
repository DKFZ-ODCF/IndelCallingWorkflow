package de.dkfz.b080.co.files;

import de.dkfz.roddy.knowledge.files.BaseFile;

/**
 * Created by heinold on 15.01.16.
 */
public class ControlBamFile extends BasicBamFile {
    public ControlBamFile(ConstructionHelperForBaseFiles helper) {
        super(helper);
    }

    public ControlBamFile(BaseFile parent) {
        super(parent);
    }

    public ControlBamFile(BasicBamFile parent) {
        super(parent);
    }
}
