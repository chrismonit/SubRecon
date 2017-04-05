package subrecon;

import com.beust.jcommander.Parameter;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class CommandArgs {
    
    // input data
    @Parameter(names = {"-sequences", "-s"}, required = true, description = "Amino acid sequence alignment (assumes FASTA format by default)")
    private String alignmentPath;
    
    public String getAlignPath(){
        return alignmentPath;
    }
    
    
    @Parameter(names = {"-tree", "-t"}, required = true, description = "Newick format tree file, rooted on the branch of interest (NB branch lengths should have been estimated in advance with the same model specified here)")
    private String treePath;
    
    public String getTreePath(){
        return treePath;
    }
    
    
    @Parameter(names = {"-phy"}, required = false, description = "Alignment is in Phylip format")
    private boolean phy = false;
    
    public boolean getPhy(){
        return phy;
    }
    
    // model 
    @Parameter(names = {"-model", "-m"}, required = true, description = "Amino acid substitution model: dayhoff, jtt, wag or blosum62")
    private String model;
        
    public String getModelID(){
        return model;
    }


    // we assume in most cases the user will want to provide their own frequencies
    @Parameter(names = {"-frequencies", "-pi"}, required = false, description = "Equilibrium frequencies for amino acids, delimited by comma in order: A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V. (Use model's original estimated frequencies by default)")
    private String piString;
    
    public double[] getFrequencies(){  
        if (piString == null) {
            return null;
        }
        
        double sum = 0.0;
        String[] piStrings = piString.split(",");
        
        if (piStrings.length != 20) {
            System.out.println("Error: if specifying amino acid frequencies, exactly 20 values must be provided and delimited by comma with no spaces");
            System.exit(1);
        }
        
        double[] pi = new double[piStrings.length];
        for (int i = 0; i < pi.length; i++) {
            pi[i] = Double.parseDouble(piStrings[i]);
            sum += pi[i];
        }
        
        if ( sum < 1.0-Constants.EPSILON || sum > 1.0+Constants.EPSILON ) {
            System.out.println("Warning: input equilibribum frequencies do not sum to 1. Normalised values will be used");
        }
        
        // normalise them anyway, because the PAL amino acid substitution models don't
        for (int i = 0; i < pi.length; i++) {
            pi[i] /= sum;
        }
        
        return pi;
    }
    
    // node labels
    @Parameter(names = {"-mother", "-I"}, required = true, description = "Label for node at the beginning of branch of interest")
    private String alphaLabel = "I";
    
    public String getAlphaLabel(){
        return alphaLabel;
    }
    
    @Parameter(names = {"-daughter", "-J"}, required = true, description = "Label for node at the end of the branch of interest")
    private String deltaLabel = "J";
    
    public String getDeltaLabel(){
        return deltaLabel;
    }
    
    // other options
    @Parameter(names = {"-site"}, required = false, description = "Single alignment column to analyse (analyse all columns by default)")
    private int site;
    
    public int getSite(){
        return site - 1; // user can give site number in non-zero based
    }

    @Parameter(names = {"-help", "-h"}, required = false, description = "Print help information and exit")
    private boolean showHelp = false;
    
    public boolean getShowHelp(){
        return showHelp;
    }
    
    @Parameter(names = {"-nosort"}, required = false, description = "Do NOT sort transition probabilities in descending order, but by amino acid canonical order instead")
    private boolean noSort = false;
    
    public boolean getNoSort(){
        return noSort;
    }
    
    @Parameter(names = {"-threshold"}, required = false, description = "Minimum probability value for a substitution to be displayed ('-threshold 0.0' will print all 400 possibilities)")
    private double threshold = 1./Math.pow(10., (double)Constants.SIG_DIGITS);
        
    public double getThreshold(){
        return threshold;
    }

    @Parameter(names = {"-verbose", "-v"}, required = false, description = "Print results for all sites, including those where Prob(I->I) >= 1 - [threshold] (these are omitted by default)")
    private boolean verbose = false;
    
    public boolean getVerbose(){
        return verbose;
    }
    
}
