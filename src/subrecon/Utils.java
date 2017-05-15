package subrecon;

import java.util.Arrays;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class Utils {
    
    private Utils(){}
    
    public static double roundDouble(double toRound, int decPlaces){
        long tenMultiple = (long)Math.pow(10.0, (double)decPlaces);
        return (double)Math.round( toRound * tenMultiple )  / tenMultiple;
    }
    
    public double[] getLnValues(double[] array){
        double[] logged = new double[array.length];
        for (int i = 0; i < array.length; i++) {
            logged[i] = Math.log(array[i]);
        }
        return logged;
    }
    
    public void lnArray(double[] array){
        for (int i = 0; i < array.length; i++) {
            array[i] = Math.log(array[i]);
        }
    }
    
    public static double getLogSumComponents(double[] mixLogLikelihoods){
        /* // Allows you to values which have been log-transformed
           mixLogLikelihoods == { ln(a), ln(b), ln(c), ln(d) }; // NB order doesn't matter, nor does number of entries. 
           ln(a+b+c+d) // where a is the largest
           = ln( a(1 + b/a + c/a + d/a) )
           = ln(a) + ln(1 + b/a + c/a + d/a)
           = ln(a) + ln( 1 + exp[ln(b)-ln(a)] + exp[ln(c)-ln(a)] + exp[ln(d)-ln(a)] )
        */
        
        if (mixLogLikelihoods.length == 1) {
            return mixLogLikelihoods[0];
        }
        
        int indexOfMaxL = -1;
        double max = Double.MIN_VALUE;
        
        for (int i = 0; i < mixLogLikelihoods.length; i++) {
            if (mixLogLikelihoods[i] > max) {
                max = mixLogLikelihoods[i];
                indexOfMaxL = i;
            }
        }
        System.out.println("indexOfMaxL "+indexOfMaxL);
        System.out.println(Arrays.toString(mixLogLikelihoods));
        double totalLogLikelihood = 0.0; // sum over mix components
        totalLogLikelihood += mixLogLikelihoods[indexOfMaxL]; // ln(a)

        double sum = 1.0;
        for (int i = 0; i < mixLogLikelihoods.length; i++) {
            if (i == indexOfMaxL) continue;
            sum += Math.exp( 
                    (mixLogLikelihoods[i]) - (mixLogLikelihoods[indexOfMaxL])
            );// exp
        }
        totalLogLikelihood += Math.log(sum);

        return totalLogLikelihood;
    }//getLogSumComponents
    
}
