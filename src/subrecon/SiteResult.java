package subrecon;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import pal.datatype.AminoAcids;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class SiteResult {
    
    private int site;
    private double marginalLnL;
    
    private List<Double> aboveThreshProbs;
    private List<String> aboveThreshSubs;
    
    private AminoAcids aa;
    
    
    
    public int getSite(){
        return site;
    }
    
    public double getMarginalLnL(){
        return marginalLnL;
    }
        
    public SiteResult(int site, double marginalLnL, double[][] branchProbs, 
            double threshold, boolean ignoreII, String delimiter, boolean sortByProb){
        
        // NB deliberately do not assign branchProbs to a field, because it is large and dont want it to persist in memory
        this.site = site;
        this.marginalLnL = marginalLnL;
        
        this.aboveThreshProbs = new ArrayList<Double>();
        this.aboveThreshSubs = new ArrayList<String>();
        
        this.aa = AminoAcids.DEFAULT_INSTANCE;
        
        // save the reisdue pairs with high probabilities
        // eg V->A, 0.99 etc
        for (int i = 0; i < branchProbs.length; i++) {
            for (int j = 0; j < branchProbs[0].length; j++) {
                
                if (ignoreII && i==j) continue;
                
                if (branchProbs[i][j] >= threshold) {
                    aboveThreshProbs.add(branchProbs[i][j]);
                    aboveThreshSubs.add(  new StringBuilder().append(aa.getChar(i)).append(aa.getChar(j)).toString()  );
                }
                
            } // for j            
        } // for i
        
        if (sortByProb) {
            sortByProb();
        }
        
    }// site result
    
    public SiteResult(int site, double marginalLnL, double[][] branchProbs, double threshold, boolean sortByProb){
        this(site, marginalLnL, branchProbs, threshold, false, Constants.DELIM, sortByProb);
    }
    
    
    private void sortByProb(){
        // bubble sort (descending)
        for (int i = 0; i < aboveThreshProbs.size() - 1; i++) {

            for (int j = 1; j < aboveThreshProbs.size() - i; j++) {
                if (aboveThreshProbs.get(j - 1) < aboveThreshProbs.get(j)) {
                    Collections.swap(aboveThreshProbs, j-1, j);
                    Collections.swap(aboveThreshSubs, j-1, j);
                }
            }// for j
        } // for i
    
    }
    
    
    @Override
    public String toString(){
        
        StringBuilder s = new StringBuilder();
        s.append(site+1); // correct for zero based
        s.append(Constants.DELIM);
        s.append( Utils.roundDouble(marginalLnL, Constants.SIG_DIGITS) );
        for (int i = 0; i < aboveThreshProbs.size(); i++) {
            s.append(Constants.DELIM);
            s.append(aboveThreshSubs.get(i));
            s.append(Constants.SUB_PROB_DELIM); 
            s.append( Utils.roundDouble(aboveThreshProbs.get(i), Constants.SIG_DIGITS));
        }

        return s.toString();
    }
    
}
