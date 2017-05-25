package subrecon.molevo;

import pal.alignment.Alignment;
import pal.alignment.SimpleAlignment;
import pal.datatype.DataTypeTool;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 * 
 * Wrapper for PAL alignment class, which adds a convenient method for getting a sequence's data
 */
public class AdvancedAlignmentAminoAcid extends SimpleAlignment {

    public AdvancedAlignmentAminoAcid(Alignment a){
        super(a);
    }
    
    public int getStateBySequenceName(String name, int site){ 
        int sequenceID = super.whichIdNumber(name);
        char stateAsChar = super.getData(sequenceID, site);
        int stateAsInt = DataTypeTool.getUniverisalAminoAcids().getState(stateAsChar); 
        
        return stateAsInt;
    }
    
    
    
}// class

