package subrecon.molevo;

import java.io.PrintWriter;
import pal.substmodel.AminoAcidModel;
import pal.util.XMLConstants;


/**
 * WAG model of amino acid evolution (S. Whelan and N. Goldman 2000)
 * This implementation uses the rate and frequency values found in the 
 * wag.dat file released with Ziheng Yang's PAML package
 * 
 * Original WAG class:
 * @author Korbinian Strimmer
 * @author Alexei Drummond
 *
 * This extension:
 * @author Christopher Monit
 */
public class WAGDotDat extends AminoAcidModel implements XMLConstants
{
    /**
     * constructor
     *
     * @param f amino acid frequencies
     */
    public WAGDotDat(double[] f)
    {
        super(f);

    }
    
    // Get numerical code describing the model type
    @Override
    public int getModelID()
    {
        return 5;
    }
    
    
    @Override
    public void report(PrintWriter out)
    {
        out.println("Model of substitution: WAG (Whelan-Goldman 2000) using values found in wag.dat file, released with Ziheng Yang's PAML distribution (ver. 4.9e)");
        out.println();
        printFrequencies(out);
    }
    
    /**
     * get the frequencies of the original data set that
     * formed the basis for the estimation of the rate matrix
     *
     * @param f array where amino acid frequencies will be stored
     */
    public static void getOriginalFrequencies(double[] f)
    { // values found in wag.dat
        f[0] = 0.0866279;
        f[1] = 0.043972;
        f[2] = 0.0390894;
        f[3] = 0.0570451;
        f[4] = 0.0193078;
        f[5] = 0.0367281;
        f[6] = 0.0580589;
        f[7] = 0.0832518;
        f[8] = 0.0244313;
        f[9] = 0.048466;
        f[10] = 0.086209;
        f[11] = 0.0620286;
        f[12] = 0.0195027;
        f[13] = 0.0384319;
        f[14] = 0.0457631;
        f[15] = 0.0695179;
        f[16] = 0.0610127;
        f[17] = 0.0143859;
        f[18] = 0.0352742;
        f[19] = 0.0708956;
    }

    /**
     * @return the frequencies of the original data set that
     * formed the basis for the estimation of the rate matrix
     */
    public static double[] getOriginalFrequencies() {
        double[] f = new double[20];
        getOriginalFrequencies(f);
        return f;
    }
    
    @Override
    public String getUniqueName() {
        return "wag.dat";
    }

    //
    // Private stuff
    //

    // WAG model of amino acid evolution
    // as found in Ziheng Yang's example wag.dat file in the PAML package
    @Override
    protected void rebuildRateMatrix(double[][] rate, double[] parameters)	{
            // Q matrix
        rate[0][1] = 0.551571;
        rate[0][2] = 0.509848;
        rate[1][2] = 0.635346;
        rate[0][3] = 0.738998;
        rate[1][3] = 0.147304;
        rate[2][3] = 5.429420;
        rate[0][4] = 1.027040;
        rate[1][4] = 0.528191;
        rate[2][4] = 0.265256;
        rate[3][4] = 0.030295;
        rate[0][5] = 0.908598;
        rate[1][5] = 3.035500;
        rate[2][5] = 1.543640;
        rate[3][5] = 0.616783;
        rate[4][5] = 0.098818;
        rate[0][6] = 1.582850;
        rate[1][6] = 0.439157;
        rate[2][6] = 0.947198;
        rate[3][6] = 6.174160;
        rate[4][6] = 0.021352;
        rate[5][6] = 5.469470;
        rate[0][7] = 1.416720;
        rate[1][7] = 0.584665;
        rate[2][7] = 1.125560;
        rate[3][7] = 0.865584;
        rate[4][7] = 0.306674;
        rate[5][7] = 0.330052;
        rate[6][7] = 0.567717;
        rate[0][8] = 0.316954;
        rate[1][8] = 2.137150;
        rate[2][8] = 3.956290;
        rate[3][8] = 0.930676;
        rate[4][8] = 0.248972;
        rate[5][8] = 4.294110;
        rate[6][8] = 0.570025;
        rate[7][8] = 0.249410;
        rate[0][9] = 0.193335;
        rate[1][9] = 0.186979;
        rate[2][9] = 0.554236;
        rate[3][9] = 0.039437;
        rate[4][9] = 0.170135;
        rate[5][9] = 0.113917;
        rate[6][9] = 0.127395;
        rate[7][9] = 0.030450;
        rate[8][9] = 0.138190;
        rate[0][10] = 0.397915;
        rate[1][10] = 0.497671;
        rate[2][10] = 0.131528;
        rate[3][10] = 0.084805;
        rate[4][10] = 0.384287;
        rate[5][10] = 0.869489;
        rate[6][10] = 0.154263;
        rate[7][10] = 0.061304;
        rate[8][10] = 0.499462;
        rate[9][10] = 3.170970;
        rate[0][11] = 0.906265;
        rate[1][11] = 5.351420;
        rate[2][11] = 3.012010;
        rate[3][11] = 0.479855;
        rate[4][11] = 0.074034;
        rate[5][11] = 3.894900;
        rate[6][11] = 2.584430;
        rate[7][11] = 0.373558;
        rate[8][11] = 0.890432;
        rate[9][11] = 0.323832;
        rate[10][11] = 0.257555;
        rate[0][12] = 0.893496;
        rate[1][12] = 0.683162;
        rate[2][12] = 0.198221;
        rate[3][12] = 0.103754;
        rate[4][12] = 0.390482;
        rate[5][12] = 1.545260;
        rate[6][12] = 0.315124;
        rate[7][12] = 0.174100;
        rate[8][12] = 0.404141;
        rate[9][12] = 4.257460;
        rate[10][12] = 4.854020;
        rate[11][12] = 0.934276;
        rate[0][13] = 0.210494;
        rate[1][13] = 0.102711;
        rate[2][13] = 0.096162;
        rate[3][13] = 0.046730;
        rate[4][13] = 0.398020;
        rate[5][13] = 0.099921;
        rate[6][13] = 0.081134;
        rate[7][13] = 0.049931;
        rate[8][13] = 0.679371;
        rate[9][13] = 1.059470;
        rate[10][13] = 2.115170;
        rate[11][13] = 0.088836;
        rate[12][13] = 1.190630;
        rate[0][14] = 1.438550;
        rate[1][14] = 0.679489;
        rate[2][14] = 0.195081;
        rate[3][14] = 0.423984;
        rate[4][14] = 0.109404;
        rate[5][14] = 0.933372;
        rate[6][14] = 0.682355;
        rate[7][14] = 0.243570;
        rate[8][14] = 0.696198;
        rate[9][14] = 0.099929;
        rate[10][14] = 0.415844;
        rate[11][14] = 0.556896;
        rate[12][14] = 0.171329;
        rate[13][14] = 0.161444;
        rate[0][15] = 3.370790;
        rate[1][15] = 1.224190;
        rate[2][15] = 3.974230;
        rate[3][15] = 1.071760;
        rate[4][15] = 1.407660;
        rate[5][15] = 1.028870;
        rate[6][15] = 0.704939;
        rate[7][15] = 1.341820;
        rate[8][15] = 0.740169;
        rate[9][15] = 0.319440;
        rate[10][15] = 0.344739;
        rate[11][15] = 0.967130;
        rate[12][15] = 0.493905;
        rate[13][15] = 0.545931;
        rate[14][15] = 1.613280;
        rate[0][16] = 2.121110;
        rate[1][16] = 0.554413;
        rate[2][16] = 2.030060;
        rate[3][16] = 0.374866;
        rate[4][16] = 0.512984;
        rate[5][16] = 0.857928;
        rate[6][16] = 0.822765;
        rate[7][16] = 0.225833;
        rate[8][16] = 0.473307;
        rate[9][16] = 1.458160;
        rate[10][16] = 0.326622;
        rate[11][16] = 1.386980;
        rate[12][16] = 1.516120;
        rate[13][16] = 0.171903;
        rate[14][16] = 0.795384;
        rate[15][16] = 4.378020;
        rate[0][17] = 0.113133;
        rate[1][17] = 1.163920;
        rate[2][17] = 0.071917;
        rate[3][17] = 0.129767;
        rate[4][17] = 0.717070;
        rate[5][17] = 0.215737;
        rate[6][17] = 0.156557;
        rate[7][17] = 0.336983;
        rate[8][17] = 0.262569;
        rate[9][17] = 0.212483;
        rate[10][17] = 0.665309;
        rate[11][17] = 0.137505;
        rate[12][17] = 0.515706;
        rate[13][17] = 1.529640;
        rate[14][17] = 0.139405;
        rate[15][17] = 0.523742;
        rate[16][17] = 0.110864;
        rate[0][18] = 0.240735;
        rate[1][18] = 0.381533;
        rate[2][18] = 1.086000;
        rate[3][18] = 0.325711;
        rate[4][18] = 0.543833;
        rate[5][18] = 0.227710;
        rate[6][18] = 0.196303;
        rate[7][18] = 0.103604;
        rate[8][18] = 3.873440;
        rate[9][18] = 0.420170;
        rate[10][18] = 0.398618;
        rate[11][18] = 0.133264;
        rate[12][18] = 0.428437;
        rate[13][18] = 6.454280;
        rate[14][18] = 0.216046;
        rate[15][18] = 0.786993;
        rate[16][18] = 0.291148;
        rate[17][18] = 2.485390;
        rate[0][19] = 2.006010;
        rate[1][19] = 0.251849;
        rate[2][19] = 0.196246;
        rate[3][19] = 0.152335;
        rate[4][19] = 1.002140;
        rate[5][19] = 0.301281;
        rate[6][19] = 0.588731;
        rate[7][19] = 0.187247;
        rate[8][19] = 0.118358;
        rate[9][19] = 7.821300;
        rate[10][19] = 1.800340;
        rate[11][19] = 0.305434;
        rate[12][19] = 2.058450;
        rate[13][19] = 0.649892;
        rate[14][19] = 0.314887;
        rate[15][19] = 0.232739;
        rate[16][19] = 1.388230;
        rate[17][19] = 0.365369;
        rate[18][19] = 0.314730;
    }
}
