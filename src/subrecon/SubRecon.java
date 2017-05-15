package subrecon;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
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
            long start = System.currentTimeMillis();
            for (int iSite = 0; iSite < alignment.getLength(); iSite++) {
                
                SiteResult result = analyseSite(iSite);
                if (verbose || result.getMaxIIProb() <= 1.-threshold) {
                    interestingSite = true; // at least one site has result to be printed
                    System.out.println(result);
                }// else print nothing
            }// for
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
        // i * pi.length + j + iRate * pi.length*pi.length
        double[][][] logConditionalMix = new double[pi.length][pi.length][gRates.getNumberOfRates()]; 
        //double sumJointStateProbs = 0.0; // will equal the total likelihood for this site, not conditional on states or rate class
        
        for (int iRate = 0; iRate < gRates.getNumberOfRates(); iRate++) {
            double rate = gRates.getRate(iRate);        
            
            // compute conditional lnls
            Count alphaScalingCount = new Count();
            Count deltaScalingCount = new Count();
            // these have not yet been corrected for scaling:
            double[] alphaConditionals = downTreeMarginal(nodeA, site, alphaScalingCount, rate); // TODO should log all of these here (and correct for scaling?)
            double[] deltaConditionals = downTreeMarginal(nodeD, site, deltaScalingCount, rate);            
            
            double logScalingCorrection = alphaScalingCount.get() + deltaScalingCount.get(); //NB this is a logged value
            //System.out.println("scalingCorrection "+scalingCorrection);
            model.setDistance(branchLengthAD * rate);

            for (int iAlpha = 0; iAlpha < pi.length; iAlpha++) {

                double logAlphaTerms = Math.log(pi[iAlpha]) + Math.log(alphaConditionals[iAlpha]); // NB could store log pi values

                for (int iDelta = 0; iDelta < pi.length; iDelta++) {
                    
                    double logScaledConditionalL = logAlphaTerms + Math.log(model.getTransitionProbability(iAlpha, iDelta)) + Math.log(deltaConditionals[iDelta]); // could log both conditionals AND transition probs in advance
                    double logConditionalL = logScaledConditionalL + logScalingCorrection; // PROBLEM. if scaled loglikelihood is high magnitude negative number, taking exp below will give zero
                    logConditionalMix[iAlpha][iDelta][iRate] = logConditionalL; // contribution from this rate class
                }// iDelta
            }// iAlpha
            
        } // for iRate
        
        // compute marginal. This could be moved into the loop above
        double[] logMarginalMix = new double[gRates.getNumberOfRates()];
        for (int iRate = 0; iRate < gRates.getNumberOfRates(); iRate++) {
            Count marginalScalingCount = new Count();
            double logScaledMarginalL = Math.log( computeTotalL(site, marginalScalingCount, gRates.getRate(iRate)) ); // not corrected for scaling
            double logMarginalL = logScaledMarginalL + marginalScalingCount.get(); // correcting for scaling
            logMarginalMix[iRate] = logMarginalL;
        }
        
        double logMarginalL = Utils.getLogSumComponents(logMarginalMix); 
        // NB this is not really the marginalLL, since we've not multiplied by invNCat
        // it is the denominator in the joint state prob equation, however
        
        double[][] jointStateProbs = new double[pi.length][pi.length];
        for (int iAlpha = 0; iAlpha < pi.length; iAlpha++) {
            for (int iDelta = 0; iDelta < pi.length; iDelta++) {
                jointStateProbs[iAlpha][iDelta] = 
                        Math.exp(Utils.getLogSumComponents(logConditionalMix[iAlpha][iDelta]) - logMarginalL                        
                        );
            }
        }
        
//        // normalise
//        double invSumBranchProbs = 1./sumJointStateProbs;
//        System.out.println("1/sum="+invSumBranchProbs);
//        
//        System.out.println("sum="+sumJointStateProbs);
//        for (int iAlpha = 0; iAlpha < pi.length; iAlpha++) {
//            for (int iDelta = 0; iDelta < pi.length; iDelta++) {
//                logJointStateProbs[iAlpha][iDelta] *= invSumBranchProbs;
//            }
//        }
                // TODO restore sanity checks
//        if (sanityCheck) {
//            // 1) check sum of conditional probs equals marginal
//            // NB the marginal and sumBranchProbs both should include multiplier 1/nCat (prob of rate class), but these cancel since we take the ratio
//            double sumMarginalL = 0.0; // sum over rate classes. Called marginal because marginalises over iAlpha and iDelta states
//            for (int iRate = 0; iRate < gRates.getNumberOfRates(); iRate++) {
//                Count marginalScalingCount = new Count();
//                double scaledMarginalLL = Math.log( computeTotalL(site, marginalScalingCount, gRates.getRate(iRate)) ); // not corrected for scaling
//                double marginalL = Math.exp( scaledMarginalLL + marginalScalingCount.get() ); // correcting for scaling
//                sumMarginalL += marginalL;            
//            }// iRate
//           
//            if (sumJointStateProbs < sumMarginalL-Constants.EPSILON || sumJointStateProbs > sumMarginalL+Constants.EPSILON)
//                throw new RuntimeException("ERROR: Sum of conditional Ls and marginal L not equal. sumMarginalL="+sumMarginalL+"; sumJointStateProbs="+sumJointStateProbs);
//            // 2) having been normalised, check probs sum to 1
//            double sum = 0.0;
//            for (int iAlpha = 0; iAlpha < pi.length; iAlpha++) {
//                for (int iDelta = 0; iDelta < pi.length; iDelta++) {
//                    sum += logJointStateProbs[iAlpha][iDelta];
//                }
//            }
//            if (sum < 1.0-Constants.EPSILON || sum > 1.0+Constants.EPSILON)
//                throw new RuntimeException("ERROR: Sum of posterior probs != 1.0. sum="+sum);
//            
//        }// sanityCheck
        
        // invNCat term cancels in when computing jointStateProbs, but must include here
        // TODO this can be more efficient:
        double siteMarginalLL = Math.log(invNCat) + logMarginalL; // marginal over alpha, delta and rate classes (ie total likelihood)
        return new SiteResult(site, siteMarginalLL, jointStateProbs, threshold, sortByProb, sigDigits);
    } // analyseSite
    
    
    

    
    
    
    
    // this is bollocks
//    private double getLogSumComponents(double[] scaledMixLogLikelihoods, double[] logScalingCorrections){
//        
//        if (scaledMixLogLikelihoods.length != logScalingCorrections.length) {
//            throw new RuntimeException("ERROR: scaledMixLogLikelihoods.length != logScalingCorrections.length");
//        }
//        
//        // bubble sort (descending)
//        for (int i = 0; i < scaledMixLogLikelihoods.length - 1; i++) {
//            for (int j = 1; j < scaledMixLogLikelihoods.length - i; j++) {
//                if (scaledMixLogLikelihoods[j - 1] < scaledMixLogLikelihoods[j]) {
//                    
//                    double tmpScaledMixLogLikelihood = scaledMixLogLikelihoods[j];
//                    scaledMixLogLikelihoods[j] = scaledMixLogLikelihoods[j-1];
//                    scaledMixLogLikelihoods[j-1] = tmpScaledMixLogLikelihood;
//
//                    double tmpLogScalingCorrections = logScalingCorrections[j];
//                    logScalingCorrections[j] = logScalingCorrections[j-1];
//                    logScalingCorrections[j-1] = tmpLogScalingCorrections;
//                }
//            }// for j
//        } // for i
//        
//        int indexOfMaxL = -1;
//        for (int i = 0; i < scaledMixLogLikelihoods.length; i++) {
//            
//        }
//        
//        double correctedLogLikelihood = 0.0;
//        correctedLogLikelihood += scaledMixLogLikelihoods[0] + logScalingCorrections[0];
//        
//        double sum = 1.0;
//        for (int i = 1; i < scaledMixLogLikelihoods.length; i++) {
//            sum += Math.exp( 
//                    (scaledMixLogLikelihoods[i] + logScalingCorrections[i]) - 
//                            (scaledMixLogLikelihoods[0] + logScalingCorrections[0])
//            );// exp
//        }
//        correctedLogLikelihood += Math.log(sum);
//                
//        return correctedLogLikelihood;
//    }
    
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
