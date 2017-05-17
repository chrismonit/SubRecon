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
    private boolean sanityCheck = true;
    private AdvancedAlignmentAminoAcid alignment;
    private Tree tree;
    
    private AminoAcidModel model;    
    private double[] pi;
    private double[] logPi;

    private int site;
    
    private Node root;
    private Node nodeA;
    private Node nodeB;
    private double branchLengthAD;
    
    private double shape; // shape parameter for gamma dist (alpha)
    private int nCat; // rate categories for gamma dist
    double logNCat;
    
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
            double totalLnL = 0.0; // across sites
            
            boolean interestingSite = false;
            long start = System.currentTimeMillis();
            for (int iSite = 0; iSite < alignment.getLength(); iSite++) {
                
                SiteResult result = analyseSite(iSite);
                totalLnL += result.getMarginalLnL();
                if (verbose || result.getMaxIIProb() <= 1.-threshold) {
                    interestingSite = true; // at least one site has result to be printed
                    System.out.println(result);
                }// else print nothing
            }// for
            System.out.println("Total lnL: "+totalLnL);

            long duration = System.currentTimeMillis() - start;
            System.out.println("duration: "+(duration/1000) + " s");
            
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
        this.logNCat = Math.log(nCat);
        

        try{ // check input parameters are ok
            loadData(comArgs.getAlignPath(), comArgs.getTreePath(), comArgs.getPhy());
            this.model = assignModel(comArgs.getModelID(), comArgs.getFrequencies());
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
        
        this.pi = this.model.getEquilibriumFrequencies();
        this.logPi = Utils.getLnValues(pi);
        
        nodeA = root.getChild(0);  
        nodeB = root.getChild(1);
        
        /* 
            The original root is along the branch connecting A and D.
            Therefore the length of the branch connecting A and D is the sum 
            of the two branches descending from the original root.
        */ 
        branchLengthAD = nodeA.getBranchLength() + nodeB.getBranchLength();
        
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
    
    
    
    
    private AminoAcidModel assignModel(String modelID, double[] frequencies) throws ParameterException {
    
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
                throw new ParameterException("ERROR: Model identifier not recognised");
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
                throw new ParameterException("ERROR: Model identifier not recognised");
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
        double[] logConditionalMix = new double[pi.length * pi.length * gRates.getNumberOfRates()]; 
        //double sumJointStateProbs = 0.0; // will equal the total likelihood for this site, not conditional on states or rate class
        
        for (int iRate = 0; iRate < gRates.getNumberOfRates(); iRate++) {
            double rate = gRates.getRate(iRate);        
            
            // compute conditional lnls
            Count alphaScalingCount = new Count();
            Count betaScalingCount = new Count();
            
            // these have not yet been corrected for scaling. TODO could we correct for scaling here?
            double[] logAlphaConditionals = Utils.getLnValues(downTreeMarginal(nodeA, site, alphaScalingCount, rate));
            double[] logBetaConditionals = Utils.getLnValues(downTreeMarginal(nodeB, site, betaScalingCount, rate));            
            
            double logScalingCorrection = alphaScalingCount.get() + betaScalingCount.get(); //NB this is a logged value
            
            model.setDistance(branchLengthAD * rate);

            for (int iAlpha = 0; iAlpha < pi.length; iAlpha++) {

                double logAlphaTerms = logPi[iAlpha] + logAlphaConditionals[iAlpha]; // NB could store log pi values

                for (int iBeta = 0; iBeta < pi.length; iBeta++) {
                    
                    double logScaledConditionalL = logAlphaTerms + Math.log(model.getTransitionProbability(iAlpha, iBeta)) + logBetaConditionals[iBeta]; // could log both conditionals AND transition probs in advance
                    double logConditionalL = logScaledConditionalL + logScalingCorrection;
                    logConditionalMix[flatIndex(iAlpha,iBeta,iRate)] = logConditionalL; // contribution from this rate class
                }// iBeta
            }// iAlpha
            
        } // for iRate
                
        double logMarginalL = Utils.getLnSumComponents(logConditionalMix); // NB this is not really the marginalLL, since we've not multiplied by 1/nCat
        
        // compute the joint probabilities we're actually interested in, by summing over rate classes
        double[][] jointStateProbs = new double[pi.length][pi.length]; // NB this certainly needs to be a new instance
        double[] logRateConditionals = new double[gRates.getNumberOfRates()];
        
        for (int iAlpha = 0; iAlpha < pi.length; iAlpha++) {
            for (int iBeta = 0; iBeta < pi.length; iBeta++) {
                
                int zone = flatIndex(iAlpha, iBeta, 0); // region of array corresponding to this combination of iAlpha and iBeta
                for (int iRate = 0; iRate < gRates.getNumberOfRates(); iRate++) {
                    logRateConditionals[iRate] = logConditionalMix[zone+iRate];
                }
                
                jointStateProbs[iAlpha][iBeta] = Math.exp(Utils.getLnSumComponents(logRateConditionals) - logMarginalL);
            }// iBeta
        }// iAlpha
        
        if (sanityCheck) {
            // 1) check sum of conditional probs equals marginal
            // NB this does not strictly compute the log marginal likelihood, since we omit the 1/nCat term, which cancels in the jointStateProbs
            double[] logComputedMarginalMix = new double[gRates.getNumberOfRates()]; // marginal likelihoods for each part of the rate mixture model
            for (int iRate = 0; iRate < gRates.getNumberOfRates(); iRate++) {
                Count marginalScalingCount = new Count();
                double logScaledMarginalL = Math.log( computeTotalL(site, marginalScalingCount, gRates.getRate(iRate)) ); // not corrected for scaling
                double logComputedMarginalL = logScaledMarginalL + marginalScalingCount.get(); // correcting for scaling
                logComputedMarginalMix[iRate] = logComputedMarginalL;
            }
            double logSumComputedMarginalL = Utils.getLnSumComponents(logComputedMarginalMix); 
           
            if (logSumComputedMarginalL < logMarginalL-Constants.EPSILON || logSumComputedMarginalL > logMarginalL+Constants.EPSILON)
                throw new RuntimeException("ERROR: Sum of conditional Ls and marginal L not equal. logSumComputedMarginalL="+logSumComputedMarginalL+"; logMarginalL="+logMarginalL);
            
            //2) having been normalised, check probs sum to 1
            double sum = 0.0;
            for (int iAlpha = 0; iAlpha < pi.length; iAlpha++) {
                for (int iBeta = 0; iBeta < pi.length; iBeta++) {
                    sum += jointStateProbs[iAlpha][iBeta];
                }
            }
            if (sum < 1.0-Constants.EPSILON || sum > 1.0+Constants.EPSILON)
                throw new RuntimeException("ERROR: Sum of posterior probs != 1.0. sum="+sum);
        }// sanityCheck
        
        // 1/nCat term cancels in when computing jointStateProbs, but must include here
        double siteMarginalLL = logMarginalL - Math.log(nCat); // marginal over alpha, beta and rate classes (ie total likelihood)
        return new SiteResult(site, siteMarginalLL, jointStateProbs, threshold, sortByProb, sigDigits);
    } // analyseSite
    
    
    private int flatIndex(int i, int j, int k){
        //return iAlpha*pi.length + iBeta + iRate*pi.length*pi.length;
        return i * pi.length * gRates.getNumberOfRates() + j * gRates.getNumberOfRates() + k;
    }
    

    // normal pruning algorithm. Used for computing marginalL in sanity check
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
