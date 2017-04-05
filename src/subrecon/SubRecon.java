package subrecon;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import java.io.FileReader;
import java.io.PrintWriter;
import pal.alignment.AlignmentReaders;
import pal.alignment.SimpleAlignment;
import pal.datatype.AminoAcids;
import pal.substmodel.AminoAcidModel;
import pal.substmodel.BLOSUM62;
import pal.substmodel.Dayhoff;
import pal.substmodel.JTT;
import pal.substmodel.WAG;
import pal.tree.Node;
import pal.tree.ReadTree;
import pal.tree.Tree;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 * 
 * 
 * 
 *    r    |------
 |----|----|d
 |         |----- 
 |
 |
 |
 |         |-----
a|---------|g
 |         |-----
 |
 |
 |
 |         |-----
 |---------|b
           |-----

 * Diagram explaining the nomenclature used for variable names
 * The root of the tree as read in from file is at node r, which is on the branch of interest,
 * connecting alpha and delta nodes. (We choose to ignore this rooting, however, 
 * except for the fact that it tells us where the branch of interest is.)
 * The conventional pruning algorithm is used to compute likelihood vectors for 
 * clades descending from beta, gamma and delta
 * 
 * 
 */
public class SubRecon {
    private AdvancedAlignmentAminoAcid alignment;
    private Tree tree;
    
    private AminoAcidModel model;    
    private double[] pi;

    private int site;
    
    private Node root;
    private Node alpha;
    private Node beta;
    private Node gamma;
    private Node delta;
    private double alphaToDeltaBL;
    
    private boolean sortByProb; // sort by value for output
    private double threshold; // minimum transition probability for printing 
    private boolean verbose;
    
    private CommandArgs comArgs;
    
    public SubRecon(){}
    
    public static void main(String[] args) {
        new SubRecon().run(args);
    }
    
    public void run(String[] args){
        
        this.init(args);
        
        if (site < 0) { // default site value is -1. User has not specified a site to analyse, so we analyse all of them
            boolean interestingSite = false;
            for (int iSite = 0; iSite < alignment.getLength(); iSite++) {
                
                SiteResult result = analyseSite(iSite);
                if (verbose || result.getMaxIIProb() <= 1.-threshold) {
                    interestingSite = true; // at least one site has result to be printed
                    System.out.println(result);
                }// else print nothing
            }// for
            
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
            System.out.println(ex.getMessage());
        }
        
        if (comArgs.getShowHelp()) {
            helpAndExit(jcom);
        }
        
        loadData(comArgs.getAlignPath(), comArgs.getTreePath(), comArgs.getPhy());
        
        this.site = comArgs.getSite(); // default value is -1, meaning analyse all sites
        this.sortByProb = !comArgs.getNoSort();
        this.threshold =  comArgs.getThreshold();
        this.verbose = comArgs.getVerbose();
        
        if (site > alignment.getLength()-1) {
            throw new RuntimeException("ERROR: -site value is greater than the number of sites in the alignment");
        }
        
        String alphaLabel = comArgs.getAlphaLabel();
        String deltaLabel = comArgs.getDeltaLabel();
        
        this.model = assignModel(comArgs.getModelID(), comArgs.getFrequencies());
        this.pi = this.model.getEquilibriumFrequencies();

        // --- initialising tree and nodes etc ---
        
        root = this.tree.getRoot();
        
        if (root.getChildCount() > 2) {
            throw new RuntimeException("ERROR: Tree root unexpectedly has more than two descendents. Is tree rooted correctly?");
        }
        
        // assign alpha and delta nodes
        alpha = delta = null;
        
        for (int iChild = 0; iChild < root.getChildCount(); iChild++) {
            if (root.getChild(iChild).getIdentifier().getName().equals(alphaLabel)) {
                alpha = root.getChild(iChild);
            } else if (root.getChild(iChild).getIdentifier().getName().equals(deltaLabel)){
                delta = root.getChild(iChild);
            }
        }
        
        if (alpha == null || delta == null) {
            throw new RuntimeException("ERROR: Failed to recognise mother and daughter nodes. Are they labelled correctly in the tree?");
        }
        
        /* In the teminology used here, alpha node is treated as the root - see diagram
            However, the original root is actually along the branch connecting alpha and delta.
            Therefore the length of the branch connecting alpha and delta is the sum 
            of the two branches descending from the original root.
        */ 
        alphaToDeltaBL = alpha.getBranchLength()+delta.getBranchLength();

        
        beta = alpha.getChild(0); // desginations of beta and gamma are arbitrary
        gamma = alpha.getChild(1);
        
        try{
            PrintWriter writer = new PrintWriter(System.out);
            model.report(writer);
            writer.flush();
        }catch(Exception e){}
        
    } // init
    
    
    
    
    
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
    
    private void helpAndExit(JCommander jcom){
        
        HelpHandler handler = new HelpHandler();
        handler.printSynopsis();
        handler.printOptions(jcom);
        System.exit(0);
    }

    
    public SiteResult analyseSite(int site){

        // compute conditional lnls
        
        Count betaScalingCount = new Count();
        Count gammaScalingCount = new Count();
        Count deltaScalingCount = new Count();
        Count marginalScalingCount = new Count();

        // these have not been corrected for scaling
        double[] betaConditionals = downTreeMarginal(beta, site, betaScalingCount);
        double[] gammaConditionals = downTreeMarginal(gamma, site, gammaScalingCount);
        double[] deltaConditionals = downTreeMarginal(delta, site, deltaScalingCount);

        // compute marginal lnl
        
        double scaledMarginalL = computeTotalL(site, marginalScalingCount); // not corrected for scaling
        double scaledMarginalLL = Math.log(scaledMarginalL);
        double marginalLL = scaledMarginalLL + marginalScalingCount.get(); // correcting for scaling
        double marginalL = Math.exp(marginalLL); // used in sanity checks below
        
        double[][] branchProbs = new double[pi.length][pi.length]; // stores the 400 transition probs for branch of interest
        
        double sumConditionalL = 0.0; // for sanity check
        
        double scalingCorrection = betaScalingCount.get() + gammaScalingCount.get() + deltaScalingCount.get();
        
        // TODO this is much more efficient for speed, at the expense of more memory use
        // can we reduce memory use somehow?
        double[][] alphaDeltaProbs = new double[pi.length][pi.length];
        double[][] alphaBetaProbs = new double[pi.length][pi.length];
        double[][] alphaGammaProbs = new double[pi.length][pi.length];
        
        model.setDistance(alphaToDeltaBL);
        model.getTransitionProbabilities(alphaDeltaProbs);
        
        model.setDistance(beta.getBranchLength());
        model.getTransitionProbabilities(alphaBetaProbs);
        
        model.setDistance(gamma.getBranchLength());
        model.getTransitionProbabilities(alphaGammaProbs);
        
        // Conditional Likelihoods
         
        for (int iAlpha = 0; iAlpha < pi.length; iAlpha++) {
            
            double alphaOnlyTerm = pi[iAlpha];
            
            for (int iDelta = 0; iDelta < pi.length; iDelta++) {
                
                double scaledConditionalL = 0.0; // conditional L of iAlpha->iDelta transition along branch connecting nodes alpha and delta
                
                 // alpha -> delta
                double alphaDeltaTerm = alphaDeltaProbs[iAlpha][iDelta] * deltaConditionals[iDelta];
                
                for (int iBeta = 0; iBeta < pi.length; iBeta++) {
                               
                    // alpha -> beta
                    double alphaBetaTerm = alphaBetaProbs[iAlpha][iBeta] * betaConditionals[iBeta];
                    
                    for (int iGamma = 0; iGamma < pi.length; iGamma++) {

                        // alpha -> gamma
                        double alphaGammaTerm = alphaGammaProbs[iAlpha][iGamma] * gammaConditionals[iGamma];
                                                
                        scaledConditionalL += alphaOnlyTerm * alphaDeltaTerm * alphaBetaTerm * alphaGammaTerm;
                        
                    }// iGamma
                }// iBeta
                
                double conditionalLL = Math.log(scaledConditionalL) + scalingCorrection;
                branchProbs[iAlpha][iDelta] = Math.exp(conditionalLL - marginalLL);
                sumConditionalL += Math.exp(conditionalLL); // for sanity checks
            }// iDelta
        }// iAlpha
        
        
        // --- sanity checks ----
        double probSum = 0.0;
        for (int i = 0; i < branchProbs.length; i++) {
            for (int j = 0; j < branchProbs[0].length; j++) {
                probSum += branchProbs[i][j];
            }
        }        
        
        if ( probSum < 1.0-Constants.EPSILON || probSum > 1.0+Constants.EPSILON ) {
            throw new RuntimeException("ERROR: Branch transition probabilities do not sum to 1.0. Instead: "+probSum);
        }
        
        if (sumConditionalL < marginalL-Constants.EPSILON || sumConditionalL > marginalL+Constants.EPSILON) {
            System.out.println("corrected marginalL "+marginalL);
            System.out.println("sum "+sumConditionalL);
            throw new RuntimeException("ERROR: Sum of conditional Ls and marginal L not equal");
        }
        
        return new SiteResult(site, marginalLL, branchProbs, threshold, sortByProb);
        
    } // analyseSites
    

    
    // normal pruning algorithm
    private double computeTotalL(int site, Count scalingCorrection){
        double sum = 0.0;
        double[] rootConditionals = downTreeMarginal( tree.getRoot(), site, scalingCorrection);

        for (int iRootState = 0; iRootState < pi.length; iRootState++) {
            sum +=  pi[iRootState] * rootConditionals[iRootState];
        }
        return sum;
    }
    
    // part of normal pruning algorithm
    private double[] downTreeMarginal(Node parent, int site, Count scalingCorrection){
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
                this.model.setDistance(child.getBranchLength());
                this.model.getTransitionProbabilities(P);

                double[] childConditionals = downTreeMarginal( parent.getChild(iChild), site, scalingCorrection);

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
    

    public void loadData(String alignmentPath, String treePath, Boolean readPhylip){
        try{
            SimpleAlignment simple;
            if (readPhylip)
                simple = new SimpleAlignment(AlignmentReaders.readPhylipClustalAlignment(new FileReader(alignmentPath), new AminoAcids()));
            else
                simple = new SimpleAlignment(AlignmentReaders.readFastaSequences(new FileReader(alignmentPath), new AminoAcids()));
            
            this.alignment = new AdvancedAlignmentAminoAcid(simple);
                                
            this.tree = new ReadTree(treePath);
        }
        catch(Exception e){
            System.out.println("ERROR: Unable to load alignment or tree file(s)");
            System.exit(1);
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
