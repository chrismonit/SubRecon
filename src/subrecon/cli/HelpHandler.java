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
package subrecon.cli;

import com.beust.jcommander.JCommander;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class HelpHandler {
    
    public HelpHandler(){
    
    }
    
    public void printSynopsis(){
        System.out.println("");
        System.out.println("Synopsis:");
        System.out.println("\"SubRecon: Ancestral Reconstruction at Adjacent Nodes in a Phylogeny\"");
        System.out.println("Please cite:  Christopher Monit and Richard A. Goldstein, manuscript in preparation");
        System.out.println("");
    }
    
    public void printOptions(JCommander jcom){
        jcom.usage();
        System.out.println("");
    }
    
    
    
}
