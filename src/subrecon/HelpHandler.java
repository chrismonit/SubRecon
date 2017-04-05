package subrecon;

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
        System.out.println("\"SubRecon: Probabilistic Reconstruction of Ancestral Substitutions\"");
        System.out.println("Please cite:  Monit, C. and Goldstein, R. A. manuscript in preparation");
        System.out.println("");
    }
    
    public void printOptions(JCommander jcom){
        jcom.usage();
        System.out.println("");
    }
    
    
    
}
