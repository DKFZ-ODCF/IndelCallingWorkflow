/*
* To change this template, choose Tools | Templates
* and open the template in the editor.
*/
package de.dkfz.b080.co.files;

import de.dkfz.roddy.knowledge.files.FileGroup;

import java.util.List;

//import pipelines.snv.SNVTools;

/**
* @author michael
*/
public class VCFFileGroupForIndels extends FileGroup {

    public VCFFileGroupForIndels(List<VCFFileForIndels> files) {
        super(files);
    }



}
