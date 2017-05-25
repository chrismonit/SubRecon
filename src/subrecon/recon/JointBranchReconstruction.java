/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package subrecon.recon;

import java.util.concurrent.Callable;
import pal.substmodel.AminoAcidModel;
import pal.substmodel.RateDistribution;
import pal.tree.Node;
import pal.tree.Tree;
import subrecon.molevo.AdvancedAlignmentAminoAcid;
import subrecon.Constants;
import subrecon.utils.Utils;

/**
 *
 * @author Christopher Monit c.monit.12@ucl.ac.uk
 */
public class JointBranchReconstruction implements Callable {
    
    private AdvancedAlignmentAminoAcid alignment;
    private Tree tree;
    private Node nodeA;
    private Node nodeB;
    private AminoAcidModel model;
    private boolean sanityCheck;
    private double[] pi; // TODO may want to only store logged pi, save memory
    private double[] logPi;
    private RateDistribution rateDist;
    private int site;
    
    private double threshold;
    private int sigDigits;
    private boolean sortByProb;
    
    public JointBranchReconstruction(AdvancedAlignmentAminoAcid alignment, Tree tree, 
                                    AminoAcidModel model,
                                    double[] pi, double[] logPi,
                                    RateDistribution rateDist,
                                    double threshold, int sigDigits, boolean sortByProb,
                                    boolean sanityCheck,
                                    int site
                                    ){
    
        this.alignment = alignment; this.tree = tree;
        this.nodeA = this.tree.getRoot().getChild(0); 
        this.nodeB = this.tree.getRoot().getChild(1);
        this.model = model;
        this.rateDist = rateDist;
        this.sanityCheck = sanityCheck;
    
        this.threshold = threshold;
        this.sigDigits = sigDigits;
        this.sortByProb = sortByProb;
        
        this.pi = pi;
        this.logPi = logPi;// TODO only need to pass in logPi (will need to convert some code below though)
        
        this.site = site;
    }
    
    @Override
    public SiteResult call(){
        double[] logConditionalMix = new double[pi.length * pi.length * rateDist.getNumberOfRates()]; 
        //double sumJointStateProbs = 0.0; // will equal the total likelihood for this site, not conditional on states or rate class
        
        for (int iRate = 0; iRate < rateDist.getNumberOfRates(); iRate++) {
            double rate = rateDist.getRate(iRate);        
            
            // compute conditional lnls
            Count alphaScalingCount = new Count();
            Count betaScalingCount = new Count();
            
            // these have not yet been corrected for scaling. TODO could we correct for scaling here?
            double[] logAlphaConditionals = Utils.getLnValues(downTreeMarginal(nodeA, site, alphaScalingCount, rate));
            double[] logBetaConditionals = Utils.getLnValues(downTreeMarginal(nodeB, site, betaScalingCount, rate));            
            
            double logScalingCorrection = alphaScalingCount.get() + betaScalingCount.get(); //NB this is a logged value
            
            model.setDistance((nodeA.getBranchLength() + nodeB.getBranchLength()) * rate);

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
        double[] logRateConditionals = new double[rateDist.getNumberOfRates()];
        
        for (int iAlpha = 0; iAlpha < pi.length; iAlpha++) {
            for (int iBeta = 0; iBeta < pi.length; iBeta++) {
                
                int zone = flatIndex(iAlpha, iBeta, 0); // region of array corresponding to this combination of iAlpha and iBeta
                for (int iRate = 0; iRate < rateDist.getNumberOfRates(); iRate++) {
                    logRateConditionals[iRate] = logConditionalMix[zone+iRate];
                }
                
                jointStateProbs[iAlpha][iBeta] = Math.exp(Utils.getLnSumComponents(logRateConditionals) - logMarginalL);
            }// iBeta
        }// iAlpha
        
        if (sanityCheck) {
            // 1) check sum of conditional probs equals marginal
            // NB this does not strictly compute the log marginal likelihood, since we omit the 1/nCat term, which cancels in the jointStateProbs
            double[] logComputedMarginalMix = new double[rateDist.getNumberOfRates()]; // marginal likelihoods for each part of the rate mixture model
            for (int iRate = 0; iRate < rateDist.getNumberOfRates(); iRate++) {
                Count marginalScalingCount = new Count();
                double logScaledMarginalL = Math.log( computeTotalL(site, marginalScalingCount, rateDist.getRate(iRate)) ); // not corrected for scaling
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
        double siteMarginalLL = logMarginalL - Math.log(rateDist.getNumberOfRates()); // marginal over alpha, beta and rate classes (ie total likelihood)
        return new SiteResult(site, siteMarginalLL, jointStateProbs, threshold, sortByProb, sigDigits);
    } // analyseSite
    
    
    private int flatIndex(int i, int j, int k){
        //return iAlpha*pi.length + iBeta + iRate*pi.length*pi.length;
        return i * pi.length * rateDist.getNumberOfRates() + j * rateDist.getNumberOfRates() + k;
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
