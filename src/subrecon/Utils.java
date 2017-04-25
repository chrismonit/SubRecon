package subrecon;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class Utils {
    
    public static double roundDouble(double toRound, int decPlaces){
        long tenMultiple = (long)Math.pow(10.0, (double)decPlaces);
        return (double)Math.round( toRound * tenMultiple )  / tenMultiple;
    }
    
}
