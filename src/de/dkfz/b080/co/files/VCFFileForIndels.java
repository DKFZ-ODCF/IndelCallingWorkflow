package de.dkfz.b080.co.files;

import de.dkfz.roddy.knowledge.files.BaseFile;
import de.dkfz.roddy.knowledge.files.Tuple2;
import de.dkfz.roddy.knowledge.methods.GenericMethod;

/**
 *
 * @author michael
 */
public class VCFFileForIndels extends BaseFile {

    /**
     * The file is created on filtering
     */
    private VCFFileForIndels somaticIndelsFile;

    private VCFFileForIndels somaticFunctionalIndelsFile;

    public VCFFileForIndels(BamFile parentFile) {
        super(parentFile);
    }

    /**
     * For i.e. somatic indel file
     * @param parentFile
     */
    public VCFFileForIndels(VCFFileForIndels parentFile) {
        super(parentFile);
    }


    public Tuple2<VCFFileForIndels, TextFile> annotate() {
        BamFile controlBam = (BamFile) this.getParentFiles().get(1);
        BamFile tumorBam = (BamFile) this.getParentFiles().get(0);
        return GenericMethod.callGenericTool(COConstants.TOOL_INDEL_ANNOTATION, tumorBam, controlBam, this);
    }

    public Tuple2<VCFFileForIndels, TextFile> deepAnnotate() {
        //TODO Handle virtual dependency to annotate(), so if annotate was called, add a dependency
        BamFile controlBam = (BamFile) this.getParentFiles().get(1);
        BamFile tumorBam = (BamFile) this.getParentFiles().get(0);
        return GenericMethod.callGenericTool(COConstants.TOOL_INDEL_DEEP_ANNOTATION, tumorBam, controlBam, this, "PIPENAME=INDEL_DEEPANNOTATION");
    }

    public void filter() {
        //TODO Handle virtual dependency to annotate() and deepannotate(), so if one was called, add a dependency the the later
        BamFile controlBam = (BamFile) this.getParentFiles().get(1);
        BamFile tumorBam = (BamFile) this.getParentFiles().get(0);
        GenericMethod.callGenericTool(COConstants.TOOL_FILTER_VCF_FILES, tumorBam, controlBam, this);
    }
}
