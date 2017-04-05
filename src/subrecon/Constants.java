package subrecon;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class Constants {
    
    public static final String NAME = "SubRecon";
    
    public static final String DAYHOFF_ID = "dayhoff";
    public static final String JTT_ID = "jtt";
    public static final String WAG_ID = "wag";
    public static final String BLOSUM62_ID = "blosum62";

    public static final double TINY_QUANTITY = 1e-10;
    public static final double EPSILON = 1e-6; // tolerance for sanity checks. Must be larger than TINY_QUANTITY
    
    public static final String DELIM = "\t";
    public static final String SUB_PROB_DELIM = ":"; // delimiter between sub codes and prob. e.g., if SUB_PROB_DELIM==":" then output is "VA:0.99"
    public static final int SIG_DIGITS = 2;
    
}
