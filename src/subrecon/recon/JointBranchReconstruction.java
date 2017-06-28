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
    private double logNCat;
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
        this.logNCat = Math.log(this.rateDist.getNumberOfRates());
        
        this.site = site;
    }
    
    // Callable interface
    @Override
    public SiteResult call(){
        return recon();
    }
    
    public SiteResult recon(){
        double[] logConditionalMix = new double[logPi.length * logPi.length * rateDist.getNumberOfRates()]; 
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

            for (int iAlpha = 0; iAlpha < logPi.length; iAlpha++) {

                double logAlphaTerms = logPi[iAlpha] + logAlphaConditionals[iAlpha];

                for (int iBeta = 0; iBeta < logPi.length; iBeta++) {
                    
                    double logScaledConditionalL = logAlphaTerms + Math.log(model.getTransitionProbability(iAlpha, iBeta)) + logBetaConditionals[iBeta];
                    double logConditionalL = logScaledConditionalL + logScalingCorrection;
                    logConditionalMix[flatIndex(iAlpha,iBeta,iRate)] = logConditionalL; // contribution from this rate class
                }// iBeta
            }// iAlpha
            
        } // for iRate
                
        double logSumConditionals = Utils.getLnSumComponents(logConditionalMix); // sum of conditional probs, ie sum over alpha, beta and rate classes (not strictly marginal L, as we've not multiplied by 1/nCat)
        
        // compute the joint probabilities we're actually interested in, by summing over rate classes
        double[][] jointStateProbs = new double[logPi.length][logPi.length]; // NB this certainly needs to be a new array instance, because it is passed to SiteResult instance
        double[] logRateConditionals = new double[rateDist.getNumberOfRates()];
        
        for (int iAlpha = 0; iAlpha < logPi.length; iAlpha++) {
            for (int iBeta = 0; iBeta < logPi.length; iBeta++) {
                
                int zone = flatIndex(iAlpha, iBeta, 0); // region of array corresponding to this combination of iAlpha and iBeta
                for (int iRate = 0; iRate < rateDist.getNumberOfRates(); iRate++) {
                    logRateConditionals[iRate] = logConditionalMix[zone+iRate];
                }
                
                jointStateProbs[iAlpha][iBeta] = Math.exp(Utils.getLnSumComponents(logRateConditionals) - logSumConditionals);
            }// iBeta
        }// iAlpha
        
        if (sanityCheck) {
            checkConditionalsSumToMarginal(logSumConditionals);
            checkSumToOne(jointStateProbs);
        }// sanityCheck
        
        double siteMarginalLL = logSumConditionals - logNCat; // marginal over alpha, beta and rate classes (ie total site likelihood). 1/nCat term cancels in when computing jointStateProbs, but must include here
        return new SiteResult(site, siteMarginalLL, jointStateProbs, threshold, sortByProb, sigDigits);
    } // analyseSite
    
    private void checkConditionalsSumToMarginal(double logSumConditionals){
        // NB this does not strictly compute the log marginal likelihood, since we omit the 1/nCat term, which cancels in the jointStateProbs
        double[] logComputedMarginalMix = new double[rateDist.getNumberOfRates()]; // marginal likelihoods for each part of the rate mixture model
        for (int iRate = 0; iRate < rateDist.getNumberOfRates(); iRate++) {
            Count marginalScalingCount = new Count();
            double logScaledMarginalL = Math.log( computeTotalL(site, marginalScalingCount, rateDist.getRate(iRate)) ); // not corrected for scaling
            logComputedMarginalMix[iRate] = logScaledMarginalL + marginalScalingCount.get(); // correcting for scaling
        }
        double logSumComputedMarginalL = Utils.getLnSumComponents(logComputedMarginalMix); 

        if (logSumComputedMarginalL < logSumConditionals-Constants.EPSILON || logSumComputedMarginalL > logSumConditionals+Constants.EPSILON)
            throw new RuntimeException("ERROR: Failed sanity check. Sum of conditional Ls and marginal L not equal. logSumComputedMarginalL="+logSumComputedMarginalL+"; logMarginalL="+logSumConditionals);
    }
        
    
    private void checkSumToOne(double[][] jointStateProbs){ // these values will not be logged
        double sum = 0.0;
        for (int iAlpha = 0; iAlpha < logPi.length; iAlpha++) {
            for (int iBeta = 0; iBeta < logPi.length; iBeta++) {
                sum += jointStateProbs[iAlpha][iBeta];
            }
        }
        if (sum < 1.0-Constants.EPSILON || sum > 1.0+Constants.EPSILON)
            throw new RuntimeException("ERROR: Failed sanity check. Sum of posterior probs != 1.0. sum="+sum);
    }
    
    
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
