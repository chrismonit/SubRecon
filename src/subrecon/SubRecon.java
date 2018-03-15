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

import subrecon.molevo.*;
import subrecon.utils.*;
import subrecon.recon.*;
import subrecon.cli.*;
import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import pal.alignment.AlignmentParseException;
import pal.alignment.AlignmentReaders;
import pal.alignment.SimpleAlignment;
import pal.datatype.AminoAcids;
import pal.substmodel.AminoAcidModel;
import pal.substmodel.BLOSUM62;
import pal.substmodel.Dayhoff;
import pal.substmodel.GammaRates;
import pal.substmodel.JTT;
import pal.substmodel.RateDistribution;
import pal.substmodel.WAG;
import pal.tree.Node;
import pal.tree.NodeUtils;
import pal.tree.ReadTree;
import pal.tree.Tree;
import pal.tree.TreeParseException;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 * 
 * Fig. 1 Representation of tree showing node nomenclature used here. r = root.
 * The conventional pruning algorithm is used to compute conditional probability 
 * vectors for clades A and B.
 * 
 *       |-----
 *   |---|A
 *   |   |-----
 *   |
 * r-|
 *   |
 *   |   |----
 *   |---|B
 *       |----
 * 
 */
public class SubRecon {
    private boolean sanityCheck;
    private AdvancedAlignmentAminoAcid alignment;
    private Tree tree;
    
    private double[] pi;
    private double[] logPi;

    private int site;
    
    private boolean sortByProb; // sort by value for output
    private double threshold; // minimum transition probability for printing 
    private boolean verbose;
    private int sigDigits;
    
    private RateDistribution rateDist;
    
    private CommandArgs comArgs;
    
    private int nThreads;
    
    public SubRecon(){}
    
    public static void main(String[] args) {
        new SubRecon().run(args);
    }
    
    public void run(String[] args){
        
        //long start = System.currentTimeMillis();
        this.init(args);
        boolean printingSites = false;
        if (site > -1) { // analysing a single site
            SiteResult result = new JointBranchReconstruction(alignment, tree, getModelInstance(comArgs.getModelID(), pi), pi, logPi, rateDist, threshold, sigDigits, sortByProb, sanityCheck, site).call();
            if (verbose || (result.getMaxIIProb() <= 1.-threshold && result.getMaxProb() >= threshold)) { // print if highest prob substitution is NOT I->I *and* if that sub has prob >= threshold
                printingSites = true; // at least one site has result to be printed
                System.out.println(result);
            }
        }else{
            ExecutorService threadPool = Executors.newFixedThreadPool(nThreads);        
            
            List<Future<SiteResult>> siteResults = new ArrayList<Future<SiteResult>>();
            for (int iSite = 0; iSite < alignment.getLength(); iSite++) {
                Future<SiteResult> siteResult = threadPool.submit(new JointBranchReconstruction(alignment, tree, getModelInstance(comArgs.getModelID(), pi), pi, logPi, rateDist, threshold, sigDigits, sortByProb, sanityCheck, iSite  )
                );// submit
                siteResults.add(siteResult);
            }// for iSite
            
            threadPool.shutdown();
            
            
            // print results
            double totalLnL = 0.0; // across sites
            for (int iSite = 0; iSite < siteResults.size(); iSite++) {
                SiteResult result = null;
                try{
                    result = siteResults.get(iSite).get();
                }catch(InterruptedException e){
                    System.out.printf("ERROR: Site %d%n", (iSite+1));
                    e.printStackTrace();
                    continue;
                }catch(ExecutionException e){
                    System.out.printf("ERROR: Site %d%n", (iSite+1));
                    e.printStackTrace();
                    continue;
                }

                totalLnL += result.getMarginalLnL();
                if (verbose || (result.getMaxIIProb() <= 1.-threshold && result.getMaxProb() >= threshold)) { // result is interesting if the highest prob substitution is NOT I->I *and* if that substitution has prob above threshold
                    printingSites = true; // at least one site has result to be printed
                    System.out.println(result);
                }
            }// for iSite
            System.out.printf("Total lnL: %.10f%n", totalLnL);

        }// else (analysing all sites)

        if (!printingSites) { // produce output if no sites are deemed interesting, to avoid confusion
            System.out.printf("0 sites have non-identical substitution probabilities greater than threshold value (threshold=%.5f)%n", threshold);
            System.out.println("The options -threshold, -nosort and -verbose can be used to control output detail");
        }// if
        //long duration = System.currentTimeMillis() - start;
        //System.out.printf("Duration: %d s (%d ms)%n", (duration/1000), duration);

    }// run
    
    
    private void init(String[] args){
        
        // assign fields
        this.comArgs = new CommandArgs();
        JCommander jcom = new JCommander(this.comArgs);
        jcom.setProgramName(Constants.NAME);
        
        try{
            jcom.parse(args);
        }
        catch(ParameterException ex){
            System.out.printf("ERROR: %s%n", ex.getMessage());
            helpAndExit(jcom, 1);
        }
        
        if (comArgs.getShowHelp()) {
            helpAndExit(jcom, 0);
        }
        
        this.site = comArgs.getSite(); // default value is -1, meaning analyse all sites
        this.sortByProb = !comArgs.getNoSort();
        this.threshold =  comArgs.getThreshold();
        this.verbose = comArgs.getVerbose();
        this.sigDigits = comArgs.getSigDigits();
        this.sanityCheck = comArgs.getDebug();
        this.nThreads = comArgs.getNThreads();
        
        double shape = comArgs.getShape();
        int nCat = comArgs.getNCat();
        
        AminoAcidModel reportModel = null;
        Node root = null;
        
        try{ // check input parameters are ok
            loadData(comArgs.getAlignPath(), comArgs.getTreePath(), comArgs.getPhy());            
            root = this.tree.getRoot();
            reportModel = getModelInstance(comArgs.getModelID(), comArgs.getFrequencies());
            
            if (sigDigits < 1 || sigDigits > 15) 
                throw new ParameterException("ERROR: -sd (significant digits) argument must be 0 < sd < 16");
            
            if (site < -1 || site > alignment.getLength()-1) // site == -1 is the default number, meaning no value has been supplied. site < -1 means the user has given a (nonsensical) negative number
                throw new ParameterException("ERROR: -site value is less than 1 or greater than the number of sites in the alignment");
           
            if (shape < 0.0)
                throw new ParameterException("ERROR: -shape value must be greater than or equal to 0.0");
            
            if (nCat < 1)
                throw new ParameterException("ERROR: -n (number of rate categories) must be 1 or higher");
            
            if (root.getChildCount() > 2) 
                throw new ParameterException("ERROR: Tree root has more than two descendents. Is the tree rooted correctly?");
            
            if (nThreads < 1) 
                throw new ParameterException("ERROR: -T (number of threads) must be 1 or higher");
        }catch (ParameterException e){
            System.out.println(e.getMessage());
            helpAndExit(jcom, 1);
        }
        
        this.pi = reportModel.getEquilibriumFrequencies();
        this.logPi = Utils.getLnValues(pi);
        
        Node nodeA = root.getChild(0);  
        Node nodeB = root.getChild(1);

        if (comArgs.getRates() != null) {
            rateDist = new CustomRates(comArgs.getRates());
        }else{
            rateDist = new GammaRates(nCat, shape);        
        }
        
        System.out.println("SubRecon: Ancestral Reconstruction of Amino Acid Substitutions Along a Branch in a Phylogeny");
        System.out.println("Please cite: C. Monit and R. A. Goldstein, submitted");
        System.out.println("");
        
        try{
            PrintWriter writer = new PrintWriter(System.out);
            reportModel.report(writer);
            rateDist.report(writer);
            writer.flush();
        }catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }

        if (sanityCheck) {
            System.out.println("");
            System.out.println("######### sanityCheck == true #########");
        }
        
        System.out.println("");
        System.out.printf("Using %d thread(s)%n", nThreads);
        System.out.println("");
        
        System.out.println("--------------------------------------------------------------");
        System.out.println("Reconstruction probabilities are presented in the form [ab:x],");
        System.out.println("meaning x is the joint probability of residue [a] being       ");
        System.out.println("present at node [A] and residue [b] being present at node [B].");        
        System.out.printf( "Node [A] has %d tips and contains taxon %s.%n", NodeUtils.getLeafCount(nodeA), getSingleTerminalNode(nodeA).getIdentifier().getName());
        System.out.printf( "Node [B] has %d tips and contains taxon %s.%n", NodeUtils.getLeafCount(nodeB), getSingleTerminalNode(nodeB).getIdentifier().getName());
        System.out.println("--------------------------------------------------------------");
        
        System.out.println(SiteResult.getHeader());
                
    } // init
    
    private Node getSingleTerminalNode(Node parent){
        Node n;
        if (parent.isLeaf()) {
            n = parent;
        }else{
            n = getSingleTerminalNode(parent.getChild(0));
        }
        return n;
    }
    
    

    private AminoAcidModel getModelInstance(String modelArgument, double[] frequencies) throws ParameterException {
    
        AminoAcidModel model;
        if (frequencies == null) { // use default model frequencies
            if (modelArgument.equals(Constants.DAYHOFF_ID)) {
                model = new Dayhoff(Dayhoff.getOriginalFrequencies());
            }else if (modelArgument.equals(Constants.JTT_ID)){
                model = new JTT(JTT.getOriginalFrequencies());
            }else if (modelArgument.equals(Constants.WAG_ID)){
                model = new WAG(WAG.getOriginalFrequencies());
            }else if (modelArgument.equals(Constants.BLOSUM62_ID)){
                model = new BLOSUM62(BLOSUM62.getOriginalFrequencies());
            }else if (modelArgument.equals(Constants.WAG_DOT_DAT)){
                model = new WAGDotDat(WAGDotDat.getOriginalFrequencies());
            }else{
                throw new ParameterException("ERROR: Model identifier not recognised");
            }            
        }else{
            if (modelArgument.equals(Constants.DAYHOFF_ID)) {
                model = new Dayhoff(frequencies);
            }else if (modelArgument.equals(Constants.JTT_ID)){
                model = new JTT(frequencies);
            }else if (modelArgument.equals(Constants.WAG_ID)){
                model = new WAG(frequencies);
            }else if (modelArgument.equals(Constants.BLOSUM62_ID)){
                model = new BLOSUM62(frequencies);
            }else if (modelArgument.equals(Constants.WAG_DOT_DAT)){
                model = new WAGDotDat(frequencies);
            }else{
                throw new ParameterException("ERROR: Model identifier not recognised");
            }             
        }
        
        return model;
    }
    
    private void helpAndExit(JCommander jcom, int exitStatus){
        
        HelpHandler handler = new HelpHandler();
        handler.printSynopsis();
        handler.printOptions(jcom);
        System.exit(exitStatus);
    }

    public void loadData(String alignmentPath, String treePath, Boolean readPhylip) throws ParameterException {
        try{
            SimpleAlignment simple;
            if (readPhylip)
                simple = new SimpleAlignment(AlignmentReaders.readPhylipClustalAlignment(new FileReader(alignmentPath), new AminoAcids()));
            else
                simple = new SimpleAlignment(AlignmentReaders.readFastaSequences(new FileReader(alignmentPath), new AminoAcids()));
            
            this.alignment = new AdvancedAlignmentAminoAcid(simple);
                                
            this.tree = new ReadTree(treePath);
        }
        catch(TreeParseException e){
            throw new ParameterException("ERROR: Unable to parse tree file: "+e.getMessage());
        }
        catch(AlignmentParseException e){
            throw new ParameterException("ERROR: Unable to parse alignment file: "+e.getMessage());
        }
        catch(FileNotFoundException e){
            throw new ParameterException("ERROR: Unable to find alignment or tree file(s): "+e.getMessage());
        }
        catch(IOException e){
            throw new ParameterException("ERROR: Unable read tree or alignment file(s): "+e.getMessage());
        }
        catch(Exception e){
            e.printStackTrace();
            throw new ParameterException("ERROR: Problem reading alignment or tree");
        }

    }//loadData
    
}
