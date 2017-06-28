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
package subrecon.molevo;

import java.io.PrintWriter;
import pal.substmodel.RateDistribution;
import subrecon.Constants;

/**
 *
 * @author Christopher Monit c.monit.12@ucl.ac.uk
 */
public class CustomRates extends RateDistribution {
    
    private String msg = "Not supported - CustomRates distribution is defined by user parameterisation only";
    
    
    public CustomRates(double[] rate, double[] probability){
        super(rate.length);
        
        if (rate.length != probability.length) {
            throw new IllegalArgumentException("Number of rates and category probabilities must be equal");
        }
        
        double z = 0.0;
        for (int i = 0; i < probability.length; i++) {
            z += probability[i];
        }
        
        if (z < 1.0-Constants.EPSILON || z > 1.0+Constants.EPSILON) {
            throw new IllegalArgumentException("Category probabilities must sum to 1.0. Sum="+z);
        }
        
        super.rate = rate;
        super.probability = probability;
    }
    
    public CustomRates(double[] rate){
        super(rate.length);

        double[] probs = new double[rate.length];
        double p = 1./(double)rate.length;
        for (int i = 0; i < probs.length; i++) {
            probs[i] = p;
        }
        super.rate = rate;
        super.probability = probs;
    }
    
    @Override
    public int getNumParameters() {
        return 0; // no free parameters
    }

    @Override
    public void setParameter(double d, int i) {
        throw new UnsupportedOperationException(msg); 
    }

    @Override
    public double getParameter(int i) {
        throw new UnsupportedOperationException(msg); 
    }

    @Override
    public void setParameterSE(double d, int i) {
        throw new UnsupportedOperationException(msg); 
    }

    @Override
    public double getLowerLimit(int i) {
        throw new UnsupportedOperationException(msg); 
    }

    @Override
    public double getUpperLimit(int i) {
        throw new UnsupportedOperationException(msg); 
    }

    @Override
    public double getDefaultValue(int i) {
        throw new UnsupportedOperationException(msg); 
    }

    @Override
    public void report(PrintWriter out)
    {
            out.println("Model of rate heterogeneity: User defined rates");
            out.println("Number of rate categories: " + numRates);
            out.println();
            printRates(out);
    }
}
