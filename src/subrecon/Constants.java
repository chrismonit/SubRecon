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
    public static final String WAG_DOT_DAT = "wag.dat";
    
    public static final double TINY_QUANTITY = 1e-10;
    public static final double EPSILON = 1e-6; // tolerance for sanity checks. Must be larger than TINY_QUANTITY
    
    public static final String DELIM = "\t";
    public static final String SUB_PROB_DELIM = ":"; // delimiter between sub codes and prob. e.g., if SUB_PROB_DELIM==":" then output is "VA:0.99"
    
    public static final double DEFAULT_PRINT_THRESHOLD = 0.4;
    public static final int DEFAULT_SIG_DIGITS = 2;
    
}
