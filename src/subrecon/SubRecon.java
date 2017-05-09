package subrecon;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import pal.alignment.AlignmentParseException;
import pal.alignment.AlignmentReaders;
import pal.alignment.SimpleAlignment;
import pal.datatype.AminoAcids;
import pal.substmodel.AminoAcidModel;
import pal.substmodel.BLOSUM62;
import pal.substmodel.Dayhoff;
import pal.substmodel.GammaRates;
import pal.substmodel.JTT;
import pal.substmodel.WAG;
import pal.tree.Node;
import pal.tree.ReadTree;
import pal.tree.Tree;
import pal.tree.TreeParseException;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 * 
 * Fig. 1 Representation of tree showing node nomenclature used here. r = root.
 * The conventional pruning algorithm is used to compute conditional probability 
 * vectors for clades A and D.
 * 
 *       |-----
 *   |---|D
 *   |   |-----
 *   |
 * r-|
 *   |
 *   |   |----
 *   |---|A
 *       |----
 * 
 */
public class SubRecon {
    private AdvancedAlignmentAminoAcid alignment;
    private Tree tree;
    
    private AminoAcidModel model;    
    private double[] pi;

    private int site;
    
    private Node root;
    private Node nodeA;
    private Node nodeD;
    private double branchLengthAD;
    
    private double shape; // shape parameter for gamma dist (alpha)
    private int nCat; // rate categories for gamma dist
    double invNCat;
    
    private boolean sortByProb; // sort by value for output
    private double threshold; // minimum transition probability for printing 
    private boolean verbose;
    private int sigDigits;
    
    private GammaRates gRates;
    
    private CommandArgs comArgs;
    
    public SubRecon(){}
    
    public static void main(String[] args) {
        new SubRecon().run(args);
    }
    
    public void run(String[] args){
        
        this.init(args);
        
        if (site < 0) { // default site value is -1. User has not specified a site to analyse, so we analyse all of them
            boolean interestingSite = false;
            //long start = System.currentTimeMillis();
            for (int iSite = 0; iSite < alignment.getLength(); iSite++) {
                
                SiteResult result = analyseSite(iSite);
                if (verbose || result.getMaxIIProb() <= 1.-threshold) {
                    interestingSite = true; // at least one site has result to be printed
                    System.out.println(result);
                }// else print nothing
            }// for
            //long duration = System.currentTimeMillis() - start;
            //System.out.println("duration: "+duration);
            
            if (!interestingSite) { // produce output if no sites are deemed interesting, to avoid confusion
                System.out.println("0 sites have non-identical substitution probabilities greater than threshold value");
                System.out.println("Threshold="+threshold);
                System.out.println("The options -threshold, -nosort and -verbose can be used to control output detail");
            }// if
        }else{
            System.out.println(analyseSite(site)); // a single named site is being analysed
        }

    }
    
    
    private void init(String[] args){

        // assign fields
        this.comArgs = new CommandArgs();
        JCommander jcom = new JCommander(this.comArgs);
        jcom.setProgramName(Constants.NAME);
        
        try{
            jcom.parse(args);
        }
        catch(ParameterException ex){
            System.out.println("ERROR: " + ex.getMessage());
            helpAndExit(jcom, 1);
        }
        
        if (comArgs.getShowHelp()) {
            helpAndExit(jcom, 0);
        }
        
        this.site = comArgs.getSite(); // default value is -1, meaning analyse all sites
        this.sortByProb = !comArgs.getNoSort();
        this.threshold =  comArgs.getThreshold();
        this.verbose = comArgs.getVerbose();
        this.sigDigits = comArgs.getSigDigits();
        
        this.shape = comArgs.getShape();
        this.nCat = comArgs.getNCat();
        this.invNCat = 1./(double)nCat;
        
        // TODO this should be moved to try block, should throw exception if model identifier makes no sense
        this.model = assignModel(comArgs.getModelID(), comArgs.getFrequencies());
        this.pi = this.model.getEquilibriumFrequencies();
        
        try{ // check input parameters are ok
            loadData(comArgs.getAlignPath(), comArgs.getTreePath(), comArgs.getPhy());
            root = this.tree.getRoot();
            
            if (sigDigits < 1 || sigDigits > 15) 
                throw new ParameterException("ERROR: -sd (significant digits) argument must be 0 < sd < 16");
            
            if (site < -1 || site > alignment.getLength()-1) // site == -1 is the default number, meaning no value has been supplied. site < -1 means the user has given a (nonsensical) negative number
                throw new ParameterException("ERROR: -site value is less than 1 or greater than the number of sites in the alignment");
           
            if (shape < 0.0)
                throw new ParameterException("ERROR: -shape value must be greater than or equal to 0.0");
            
            if (nCat < 1)
                throw new ParameterException("ERROR: -n (number of rate categories) must be 1 or higher");
            
            if (root.getChildCount() > 2) 
                throw new ParameterException("ERROR: Tree root has more than two descendents. Is the tree rooted correctly?");
        
        }catch (ParameterException e){
            System.out.println(e.getMessage());
            helpAndExit(jcom, 1);
        }
        
        nodeA = root.getChild(0);  
        nodeD = root.getChild(1);
        
        /* 
            The original root is along the branch connecting A and D.
            Therefore the length of the branch connecting A and D is the sum 
            of the two branches descending from the original root.
        */ 
        branchLengthAD = nodeA.getBranchLength() + nodeD.getBranchLength();
        
        gRates = new GammaRates(nCat, shape);        
        
        try{
            PrintWriter writer = new PrintWriter(System.out);
            model.report(writer);
            gRates.report(writer);
            writer.flush();
        }catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        System.out.println("");
        
    } // init
    
    
    
    
    // TODO should throw parameter exception
    private AminoAcidModel assignModel(String modelID, double[] frequencies){
    
        AminoAcidModel model;
        if (frequencies == null) { // use default model frequencies
            if (modelID.equals(Constants.DAYHOFF_ID)) {
                model = new Dayhoff(Dayhoff.getOriginalFrequencies());
            }else if (modelID.equals(Constants.JTT_ID)){
                model = new JTT(JTT.getOriginalFrequencies());
            }else if (modelID.equals(Constants.WAG_ID)){
                model = new WAG(WAG.getOriginalFrequencies());
            }else if (modelID.equals(Constants.BLOSUM62_ID)){
                model = new BLOSUM62(BLOSUM62.getOriginalFrequencies());
            }else{
                throw new RuntimeException("ERROR: Model identifier not recognised");
            }            
        }else{
            if (modelID.equals(Constants.DAYHOFF_ID)) {
                model = new Dayhoff(frequencies);
            }else if (modelID.equals(Constants.JTT_ID)){
                model = new JTT(frequencies);
            }else if (modelID.equals(Constants.WAG_ID)){
                model = new WAG(frequencies);
            }else if (modelID.equals(Constants.BLOSUM62_ID)){
                model = new BLOSUM62(frequencies);
            }else{
                throw new RuntimeException("ERROR: Model identifier not recognised");
            }             
        }
        
        return model;
    }
    
    private void helpAndExit(JCommander jcom, int exitStatus){
        
        HelpHandler handler = new HelpHandler();
        handler.printSynopsis();
        handler.printOptions(jcom);
        System.exit(exitStatus);
    }

    
    public SiteResult analyseSite(int site){
        
        double[][] branchProbs = new double[pi.length][pi.length]; // stores the 400 posterior probs for branch of interest
        double siteMarginalL = 0.0; // this is summing over rate classes
        
        for (int iRate = 0; iRate < gRates.getNumberOfRates(); iRate++) {
            double rate = gRates.getRate(iRate);        
            
            // compute conditional lnls

            Count alphaScalingCount = new Count();
            Count deltaScalingCount = new Count();

            // these have not been corrected for scaling
            double[] alphaConditionals = downTreeMarginal(nodeA, site, alphaScalingCount, rate);
            double[] deltaConditionals = downTreeMarginal(nodeD, site, deltaScalingCount, rate);

            // compute marginal lnl (conditional on rate but not conditional on alpha or delta)
            Count marginalScalingCount = new Count();
            double scaledMarginalL = computeTotalL(site, marginalScalingCount, gRates.getRate(iRate)); // not corrected for scaling
            double scaledMarginalLL = Math.log(scaledMarginalL);

            double marginalLL = scaledMarginalLL + marginalScalingCount.get(); // correcting for scaling
            double marginalL = Math.exp(marginalLL); // used in sanity checks below
            siteMarginalL += marginalL;
            
            double sumConditionalL = 0.0; // for sanity check. Should sum to marginalL

            double scalingCorrection = alphaScalingCount.get() + deltaScalingCount.get();

            model.setDistance(branchLengthAD * rate);

            double sumBranchProb = 0.0; // for sanity check. Sum over alpha and delta for this rate class should be 1.0
            for (int iAlpha = 0; iAlpha < pi.length; iAlpha++) {

                double alphaTerms = pi[iAlpha] * alphaConditionals[iAlpha];

                for (int iDelta = 0; iDelta < pi.length; iDelta++) {

                    double scaledConditionalL = alphaTerms * model.getTransitionProbability(iAlpha, iDelta) * deltaConditionals[iDelta];

                    double conditionalLL = Math.log(scaledConditionalL) + scalingCorrection;
                    double branchProb = Math.exp(conditionalLL - marginalLL); // for this given rate class
                    branchProbs[iAlpha][iDelta] += branchProb; // sum for each of the rate classes
                    sumBranchProb += branchProb; // for sanity check
                    sumConditionalL += Math.exp(conditionalLL); // for sanity checks
                }// iDelta
            }// iAlpha

            // --- sanity checks ----
            // Sum over alpha and delta for this rate class should be 1.0
            if ( sumBranchProb < 1.0-Constants.EPSILON || sumBranchProb > 1.0+Constants.EPSILON )
                throw new RuntimeException("ERROR: Branch transition probabilities do not sum to 1.0. Instead: "+sumBranchProb+"; rateclass="+iRate);
            
            if (sumConditionalL < marginalL-Constants.EPSILON || sumConditionalL > marginalL+Constants.EPSILON)
                throw new RuntimeException("ERROR: Sum of conditional Ls and marginal L not equal. Corrected marginalL="+marginalL+"; sum="+sumConditionalL+"; rateclass="+iRate);
            // -----------------
            
            
            // because we marginalise over rate classes, we need to divide by prob of each class. These are all 1/nCat (uniform prob)
            for (int iAlpha = 0; iAlpha < 10; iAlpha++) {
                for (int iDelta = 0; iDelta < 10; iDelta++) {
                    branchProbs[iAlpha][iDelta] *= invNCat;
                }
            }
            
        } // for iRate
        
        double siteMarginalLL = Math.log(invNCat * siteMarginalL); // marginal over alpha, delta and rate classes
            
        return new SiteResult(site, siteMarginalLL, branchProbs, threshold, sortByProb, sigDigits);
        
    } // analyseSite
    

    
    // normal pruning algorithm
    private double computeTotalL(int site, Count scalingCorrection, double rate){
        double sum = 0.0;
        double[] rootConditionals = downTreeMarginal( tree.getRoot(), site, scalingCorrection, rate);

        for (int iRootState = 0; iRootState < pi.length; iRootState++) {
            sum +=  pi[iRootState] * rootConditionals[iRootState];
        }
        return sum;
    }
    
    // part of normal pruning algorithm
    private double[] downTreeMarginal(Node parent, int site, Count scalingCorrection, double rate){
        double[] parentConditionals = new double[pi.length];

        if (parent.isLeaf()){ // 'parent' is terminal node, i.e. has no children. Code here is no different from normal pruning algorithm

            String parentName = parent.getIdentifier().getName();
            int state = alignment.getStateBySequenceName(parentName, site);
            
            if (state >= 0 && state < pi.length){ //the observed state is recognised as an amino acid
                parentConditionals[state] = 1.0;
            }else  {
                for (int i = 0; i < parentConditionals.length; i++) {
                    parentConditionals[i] = 1.0; //observed state is not recognised as amino acid (may be gap). Treated as missing data. All conditional probabilities = 1.0.
                }
            }
        } else{ // NOT LEAF
    
            for (int i = 0; i < parentConditionals.length; i++) {
                parentConditionals[i] = 1.0; // multiplicative identity
            }

            for (int iChild = 0; iChild < parent.getChildCount(); iChild++){
                Node child = parent.getChild(iChild);

                double[][] P = new double[pi.length][pi.length];
                this.model.setDistance(child.getBranchLength() * rate);
                this.model.getTransitionProbabilities(P);

                double[] childConditionals = downTreeMarginal( parent.getChild(iChild), site, scalingCorrection, rate);

                for (int iParentState = 0; iParentState < pi.length; iParentState++){ // same as normal pruning algorithm
                    double sum = 0.0; //prob of observing data below this node, if the state at this node were iParentState
                    for (int jChildState = 0; jChildState < pi.length; jChildState++){
                        sum += P[iParentState][jChildState] * childConditionals[jChildState];

                    }
                    parentConditionals[iParentState] *= sum;
                } // for iParentState

            } // for iChild

        }// else (is not a leaf
        
        // scaling conditional likelihoods to prevent underflow errors
        // find biggest value
        double biggestValue = 0.0;
        for (int iParentState = 0; iParentState < pi.length; iParentState++) {
                biggestValue = Math.max(biggestValue, parentConditionals[iParentState]);
        }
        // express conditionals relative to largest value
        for (int iParentState = 0; iParentState < pi.length; iParentState++) {
                parentConditionals[iParentState] /= (biggestValue + Constants.TINY_QUANTITY); // The added small value is to prevent divide by zero problems
        }
        scalingCorrection.add( Math.log(biggestValue + Constants.TINY_QUANTITY) ); // keep track of the scaling amount as you go
        
        return parentConditionals;
    }// downTreeConditional
    

    public void loadData(String alignmentPath, String treePath, Boolean readPhylip) throws ParameterException {
        try{
            SimpleAlignment simple;
            if (readPhylip)
                simple = new SimpleAlignment(AlignmentReaders.readPhylipClustalAlignment(new FileReader(alignmentPath), new AminoAcids()));
            else
                simple = new SimpleAlignment(AlignmentReaders.readFastaSequences(new FileReader(alignmentPath), new AminoAcids()));
            
            this.alignment = new AdvancedAlignmentAminoAcid(simple);
                                
            this.tree = new ReadTree(treePath);
        }
        catch(TreeParseException e){
            throw new ParameterException("ERROR: Unable to parse tree file: "+e.getMessage());
        }
        catch(AlignmentParseException e){
            throw new ParameterException("ERROR: Unable to parse alignment file: "+e.getMessage());
        }
        catch(FileNotFoundException e){
            throw new ParameterException("ERROR: Unable to find alignment or tree file(s): "+e.getMessage());
        }
        catch(IOException e){
            throw new ParameterException("ERROR: Unable read tree or alignment file(s): "+e.getMessage());
        }
        catch(Exception e){
            e.printStackTrace();
            throw new ParameterException("ERROR: Problem reading alignment or tree");
        }

    }//loadData
    
    /*
        Keeps track of the scalling factor used in the pruning algorithm
    */
    private class Count{
        
        private double count = 0.0;
        
        public void add(double addition){
            this.count += addition;
        }
        
        public double get(){
            return this.count;
        }
        
    }// Count class
    
}
