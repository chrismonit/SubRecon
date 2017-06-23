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
