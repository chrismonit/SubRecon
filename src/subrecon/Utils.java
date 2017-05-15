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
    
    public static double[] getLnValues(double[] array){
        double[] logged = new double[array.length];
        for (int i = 0; i < array.length; i++) {
            logged[i] = Math.log(array[i]);
        }
        return logged;
    }
    
    public static void lnArray(double[] array){
        for (int i = 0; i < array.length; i++) {
            array[i] = Math.log(array[i]);
        }
    }
        
    public static double getLnSumComponents(double[] lnValues){
        /* // Allows you to values which have been log-transformed
           mixLogLikelihoods == { ln(a), ln(b), ln(c), ln(d) }; // NB order doesn't matter, nor does number of entries. 
           ln(a+b+c+d) // where a is the largest
           = ln( a(1 + b/a + c/a + d/a) )
           = ln(a) + ln(1 + b/a + c/a + d/a)
           = ln(a) + ln( 1 + exp[ln(b)-ln(a)] + exp[ln(c)-ln(a)] + exp[ln(d)-ln(a)] )
        */
        
        if (lnValues.length == 1) {
            return lnValues[0];
        }
        
        int indexOfMaxL = -1;
        double max = Double.MIN_VALUE;
        
        for (int i = 0; i < lnValues.length; i++) {
            if (lnValues[i] > max) {
                max = lnValues[i];
                indexOfMaxL = i;
            }
        }
        
        double lnSumValues = lnValues[indexOfMaxL]; // ln(a)

        double sum = 1.0;
        for (int i = 0; i < lnValues.length; i++) {
            if (i == indexOfMaxL) continue;
            sum += Math.exp((lnValues[i]) - (lnValues[indexOfMaxL]));
        }
        lnSumValues += Math.log(sum);

        return lnSumValues;
    }//getLogSumComponents
    
    // for testing
    public static void main(String[] args){
    
        double[] values = new double[]{  1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0  };
        //double[] values = new double[]{  1e-20  };
        
        for (int i = 0; i < values.length; i++) {
            values[i] += 1e-15;
        }
        
        double[] logged = Utils.getLnValues(values);
        
        System.out.println("values: "+Arrays.toString(values));
        System.out.println("logged: "+Arrays.toString(logged));
        
        double sumValues = 0.0;
        for (int i = 0; i < values.length; i++) {
            sumValues += values[i];
        }
        System.out.println("ln(sum): "+Math.log(sumValues));
        System.out.println("result : "+getLnSumComponents(logged));
                
    }
    
}
