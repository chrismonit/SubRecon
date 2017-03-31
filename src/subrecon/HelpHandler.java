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
        System.out.println("synopsis");
        System.out.println("");
    }
    
    public void printOptions(JCommander jcom){
        jcom.usage();
        System.out.println("");
    }
    
    
    
}
