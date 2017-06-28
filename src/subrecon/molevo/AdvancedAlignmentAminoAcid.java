/*
   Copyright 2017 Christopher Monit

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and 
   limitations under the License.
*/
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

