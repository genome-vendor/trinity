
#include <string>
#include <unistd.h>
#include <stdio.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <omp.h>

#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "analysis/TranscriptomeGraph.h"
#include "analysis/DNAVector.h"
#include "analysis/CompMgr.h"


const int KMER_SIZE = 24;
typedef vector<string> string_vec;
static bool NO_CLEANUP = false;


void Execute(const char * command) {
    int ret = system(command);
    if (ret != 0) {
        cout << "COMMAND: " << command << endl;
        cout << "Died with exit code " << ret << endl;
        cout << "Exiting." << endl;
        exit(-1);
    }
    
}

bool Exists(const string & s) 
{
    FILE * p = fopen(s.c_str(), "r");
    if (p != NULL) {
        fclose(p);
        return true;
    }
    // cout << "FATAL ERROR: Could not open file for read: " << s << endl;
    // cout << "Please make sure to enter the correct file name(s). Exiting now." << endl;
    
    return false;
}



void write_deBruijn_graphs(vector<string_vec>& bundled, vector<int>& component_ids, ComponentFileMgr& mgr, string outDir, string execDir, bool sStrand, map<int,bool>& pursue_component) {
    

    string cpp_graph_writer = execDir + "/../Inchworm/bin/FastaToDeBruijn";

    int num_finished = 0;
    int num_components = bundled.size();

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < bundled.size(); i++) {

        int component_id = component_ids[i];
        
        if (pursue_component.find(component_id) == pursue_component.end()) {
            // not pursuing it.
            continue;
        }


        string_vec bundled_iworms = bundled[i];
        
        string component_file_basename = mgr.GetFileName(component_id, "");
        
        
        string graph_filename = component_file_basename + ".raw.graph";
        string comp_bundle_file = component_file_basename + ".iworm_bundle";
        
        // use reads instead if iworm bundles for building de bruijn graph:
        string reads_file = component_file_basename + ".raw.fasta";
        
        stringstream cpp_graph_cmd;
        cpp_graph_cmd << cpp_graph_writer << "  --fasta " << comp_bundle_file 
                   // << "," << reads_file // use both reads and iworm contigs for building graph.
                      << " -K " << KMER_SIZE << " -C " << component_id;
        if (sStrand) {
           cpp_graph_cmd << " --SS ";
        }
        cpp_graph_cmd << " > " << graph_filename;
        

        Execute(cpp_graph_cmd.str().c_str());
        
        // cleanup
        if (! NO_CLEANUP) {
            // bundle file served its purpose.
            string cleanup_cmd = "rm " + comp_bundle_file;
            Execute(cleanup_cmd.c_str());
        }
        
        
        #pragma omp atomic
        num_finished++;
        
        #pragma omp critical 
        {

            cerr << "\r" << num_finished << "/" << num_components << " = " << (float)num_finished/num_components*100 << "% of de Bruijn graphs constructed   ";
            //cerr << "Executing: " << graph_cmd.str() << endl;
        }

    }
    
    
    cerr << "-done writing de Bruijn graphs." << endl;
    
    return;
}



void write_iworm_bundle (string filename, vector<string_vec>& bundled, vector<int>& component_ids, ComponentFileMgr& mgr) {

    ofstream ofh;
    ofh.open(filename.c_str());
    
    ofstream bundle_listing_ofh;
    bundle_listing_ofh.open("iworm_bundle_file_listing.txt");


    for (int i = 0; i < bundled.size(); i++) {

        int component_id = component_ids[i];
        string_vec bundled_iworms = bundled[i];

        stringstream s;

        s << ">s_" << component_id << endl;
        
        int num_iworms = bundled_iworms.size();
        
        for (int j = 0; j < num_iworms; j++) {
            s << bundled_iworms[j];
            if (j < num_iworms-1) {
                s << "X";
            }
        }
        s << endl;
        
        ofh << s.str();
        
        // write a local version of it according to the component identifier
        string comp_bundle_file = mgr.GetFileName(component_id, ".iworm_bundle");
        ofstream bundle_ofh;
        bundle_ofh.open(comp_bundle_file.c_str());
        bundle_ofh << s.str();
        bundle_ofh.close();

        bundle_listing_ofh << comp_bundle_file.c_str() << endl;

        
    }

    ofh.close();

    bundle_listing_ofh.close();
    
    return;
    
    
}

string get_seq_string(DNAVector& d) {

    stringstream s;
    
    for (int i = 0; i < d.isize(); i++) {
        s << d[i];
    }

    return(s.str());
}





int main(int argc,char** argv)
{
    char execPath[512];
    strcpy(execPath, argv[0]);
    execPath[strlen(execPath)-9] = 0;
    string exec = execPath;
    cout << "Path to executables: " << exec << endl;
    
    if (exec == "")
        exec = "./";
    commandArg<string> iStringCmmd("-i","read fasta file");
    commandArg<string> iwormStringCmmd("-iworm","inchworm file", "");
    commandArg<string> oStringCmmd("-o","output directory");
    commandArg<bool> pairsStringCmmd("-paired", "paired-end reads are used.", false);
    commandArg<string> butterflyCmmd("-butterfly","butterfly executable", "../Butterfly/Butterfly.jar");
    commandArg<bool> skipCmmd("-skip","skip initial 2 steps", false);
    commandArg<bool> strandCmmd("-strand","strand-specific data", false);
    commandArg<bool> nobreakCmmd("-nobreak","skip breaking", false);
    commandArg<int> minCmmd("-min","minimum sequence length", 300);
    commandArg<int> cpuCmmd("-cpu","number of CPUs to use", 10);
    commandArg<int> distCmmd("-dist","size of a read pair insert", 350);
    commandArg<int> minDumpLenCmmd("-min_all","skip components for which all seqs < min_all", 110);
    commandArg<bool> buttCmmd("-run_butterfly","runs butterfly locally", false);
    commandArg<int>  weldmerSizeCmmd("-weldmer_size", "size of weldmer to be used by GraphFromFasta", 48);
    commandArg<long> maxReadsCmd("-max_reads", "max number of reads to map to each graph", -1);
    commandArg<long> mmrCmmd("-max_mem_reads","Maximum number of reads to load into memory", -1);
    commandArg<int>  pctReadMap("-min_pct_read_mapping", "minimum percent of a read's kmers that must map to an inchworm bundle", 0);
    commandArg<int> minGlueCmmd("-min_glue", "minimum read support for glue in GraphToFasta", 2);
    commandArg<double> glueFactorCmmd("-glue_factor", "fraction of max (iworm pair coverage) for read glue support", 0.05);
    commandArg<double> minIsoRatioCmmd("-min_iso_ratio", "min ratio of (iworm pair coverage) for join", 0.05);
    commandArg<bool> debugCmmd("-debug", "verbosely describes operations", false);
    commandArg<bool> no_cleanupCmmd ("-no_cleanup", "retain input files on success", false);
    
    commandLineParser P(argc,argv);
    P.SetDescription("Assemble transcriptomes from reads.");
    P.registerArg(iStringCmmd);
    P.registerArg(iwormStringCmmd);
    P.registerArg(butterflyCmmd);
    P.registerArg(oStringCmmd);
    P.registerArg(pairsStringCmmd);
    P.registerArg(skipCmmd);
    P.registerArg(strandCmmd);
    P.registerArg(nobreakCmmd);
    P.registerArg(minCmmd);
    P.registerArg(cpuCmmd);
    P.registerArg(buttCmmd);
    P.registerArg(distCmmd);
    P.registerArg(minDumpLenCmmd);
    P.registerArg(weldmerSizeCmmd);
    P.registerArg(maxReadsCmd);
    P.registerArg(mmrCmmd);
    P.registerArg(pctReadMap);
    P.registerArg(minGlueCmmd);
    P.registerArg(glueFactorCmmd);
    P.registerArg(minIsoRatioCmmd);
    P.registerArg(debugCmmd);
    P.registerArg(no_cleanupCmmd);
    

    P.parse();

    cerr << "-------------------" << endl 
         << "---- Chrysalis ----" << endl
         << "-------------------" << endl << endl;
    

    string readString = P.GetStringValueFor(iStringCmmd);
    string iwormString = P.GetStringValueFor(iwormStringCmmd);
    string outDir = P.GetStringValueFor(oStringCmmd);
    string butterflyExec = P.GetStringValueFor(butterflyCmmd);
    

    bool PAIRED_READS_MODE = P.GetBoolValueFor(pairsStringCmmd);
    bool bSkip = P.GetBoolValueFor(skipCmmd);
    bool sStrand = P.GetBoolValueFor(strandCmmd);
    int minDumpLen = P.GetIntValueFor(minDumpLenCmmd);
    int minLen = P.GetIntValueFor(minCmmd);
    int nCPU = P.GetIntValueFor(cpuCmmd);
    int pairDist = P.GetIntValueFor(distCmmd);
   
    NO_CLEANUP = P.GetBoolValueFor(no_cleanupCmmd);
    
    bool bBreak = true;
    if (P.GetBoolValueFor(nobreakCmmd))
        bBreak = false;
    
    double glue_factor = P.GetDoubleValueFor(glueFactorCmmd);
    bool bButt = P.GetBoolValueFor(buttCmmd);
    
    long max_reads = P.GetLongValueFor(maxReadsCmd);    
    long max_mem_reads = P.GetLongValueFor(mmrCmmd);
    
    int min_pct_read_mapping = P.GetIntValueFor(pctReadMap);
    int min_glue = P.GetIntValueFor(minGlueCmmd);
    double min_iso_ratio = P.GetDoubleValueFor(minIsoRatioCmmd);
    int weldmer_size = P.GetIntValueFor(weldmerSizeCmmd);
    

    bool DEBUG = P.GetBoolValueFor(debugCmmd);
    

    if (minDumpLen > minLen) {
        minDumpLen = minLen;
    }
    
    string command;
    
    command = "mkdir ";
    command += outDir;
    system(command.c_str());
    

    /////////////////////////////////////////////////////////////////////////////////////////////////
    // Scaffolding of Inchworm contigs:
    // -run BWA to align reads to inchworm contigs, and determine scaffolding support based on pairs
    if (PAIRED_READS_MODE) {
        // Run bowtie to map reads to inchworm transcripts

        string scaffolds_finished_file = "iworm_scaffolds.txt.finished";
        if (Exists(scaffolds_finished_file) && Exists("iworm_scaffolds.txt")) {
            cerr << "Warning, Scaffolds file: iworm_scaffolds.txt previously generated and being reused." << endl;
        }
        else {
                    
            stringstream cmdstr;
            cmdstr << exec 
                   << "/../util/alignReads.pl "
                   << " --target " << iwormString
                   << " --aligner bowtie "
                   << " --seqType fa "
                   << " --single " << readString
                   << " -o iworm_bowtie "
                   << " --retain_SAM_files ";
            
            string bowtie_sam_file = "iworm_bowtie/iworm_bowtie.nameSorted.sam";
            
            if (sStrand) {
                cmdstr << " --SS_lib_type F ";
                bowtie_sam_file = "iworm_bowtie/iworm_bowtie.nameSorted.sam.+.sam";
            }
            
            // add custom options to aligner
            cmdstr << " -- -a -m 20 --best --strata --threads " << nCPU << " --quiet ";
            cmdstr << " --chunkmbs 512 ";
            
            cerr << "CMD: " << cmdstr.str() << endl;
            Execute(cmdstr.str().c_str());
            
            // Generate inchworm pair links:
            cmdstr.str(""); // clear it out.
            
            cmdstr << exec << "/../util/scaffold_iworm_contigs.pl " << bowtie_sam_file << " " << iwormString << " > iworm_scaffolds.txt";
            
            
            cerr << "CMD: " << cmdstr.str() << endl;
            
            Execute(cmdstr.str().c_str());

            string complete_scaffolds_cmd = "touch " + scaffolds_finished_file;
            Execute(complete_scaffolds_cmd.c_str());
            

            // cleanup
            if (! NO_CLEANUP) {
                // purge the alignments
                string rm_cmd = "rm -rf iworm_bowtie/ ";
                Execute(rm_cmd.c_str());
            }
            
        }
    }
    

    /////////////////////////////////////////////////////////////////////////////////////
    // GraphFromFasta:
    // Pool sequences (clustering of Inchworm contigs, join using read-support)
    // Create output file: GraphFromFasta.out
    
    cerr << "-- running Chryalis: GraphFromFasta --" << endl;


    string graphFromFastaCompleteFile = outDir + "/GraphFromIwormFasta.finished";
    
    stringstream cmdstr;

    cmdstr << exec << "GraphFromFasta -i " << iwormString
           << " -r " << readString
           << " -min_glue " <<  min_glue
           << " -glue_factor " << glue_factor
           << " -min_iso_ratio " <<  min_iso_ratio
           << " -kk " << weldmer_size;
    
    

    if (sStrand) {
        cmdstr << " -strand ";
    }
    
    cmdstr << " -report_welds ";
    //cmdstr << " -no_welds "; // disabling for now

    if (PAIRED_READS_MODE) {
        
        cmdstr << " -scaffolding iworm_scaffolds.txt ";

    }
    

    cmdstr << " > ";
    string components_file = outDir + "/GraphFromIwormFasta.out";
    cmdstr << components_file;
    
    command = cmdstr.str();
    
    if (Exists(components_file) && Exists(graphFromFastaCompleteFile)) {
        cerr << "File: " << components_file << " already exists. Using existing file assuming resume mode" << endl << endl;
    }
    else {
        cout << "Running: " << command << endl;
        Execute(command.c_str());

        // make completion-marking file
        string complete_file_cmd = "touch " + graphFromFastaCompleteFile;
        Execute(complete_file_cmd.c_str());
        
    }
    //
    ///////////////////////////////////////////////////////////////////////////////////////
    
    

    ///////////////////////////////////////////////////////////////////////////////////////
    // Read components.out, create chrysalis/bundled.fasta which represents clustered inchworm sequences
    // Create de Bruijn graphs based on clustered Inchworm seqs (partitioned chrysalis/RawComps.\d+/comp\d+.raw.graph files) 

    cerr << "-- writing inchworm bundled.fasta and computing & partitioning component graph files for parallel processing." << endl;
    
    //---- Note: This section of code requires that you have the stacksize set to unlimited. --------//

    string iString = components_file;
    //string oString = outDir + "/graph.out";
    
    string tmpName = outDir + "/tmp.fasta";
    //FILE * p = fopen(oString.c_str(), "w");
    //fclose(p);
    
    FlatFileParser parser;  
    parser.Open(iString);  // read components.out file
   
    int component_no = 0;

    ComponentFileMgr mgr(outDir);
        
    vecDNAVector tmpSeq;
    vector<int> component_ids;
    
    vector<string_vec> bundled;
    //bundled.reserve(1000000);
    string bundledName = outDir + "/bundled.fasta";
    //DNAVector separator;
    //separator.SetFromBases("X");
    string separator = "X";
    
    FILE * pOut = NULL;
        
    while (parser.ParseLine()) {
        if (parser.GetItemCount() == 0)
            continue;
        
        if (parser.AsString(0)[0] == '#') {
            continue; // ignore comments
        }
                
        if (parser.AsString(0) == "COMPONENT") {
            int component_no = parser.AsInt(1);
            int num_iworm_contigs = parser.AsInt(2);
           
            tmpSeq.resize(0); // hold the iworm seqs corresponding to a single component
            
            while (parser.ParseLine()) {
                
                if (parser.GetItemCount() == 0)
                    continue;
                if (parser.AsString(0)[0] == '#') {
                    continue; // ignore comments
                }         
                
       
                if (parser.AsString(0) == "END") {
                    break;	  
                }
                
                const char * ppp = parser.AsString(0).c_str();
                
                if (ppp[0] == '>') {
                    DNAVector tmp;
                    tmpSeq.push_back(tmp);
                    continue;
                }
                
                DNAVector app;
                app.SetFromBases(parser.AsString(0));
                //cerr << "adding: " << parser.AsString(0) << endl;
                tmpSeq[tmpSeq.isize()-1] += app; // adding the iworm sequence to that iworm entry
                
                //fprintf(p, "%s\n", parser.Line().c_str());
            }
            //fclose(p);
            
            if (tmpSeq.isize() == 1 && tmpSeq[0].isize() < minLen) {
              //cerr << "-discarding entry, too short." << endl;
                continue;
            }
            
            bool bGood = false;
            for (int x=0; x<tmpSeq.isize(); x++) {
                if (tmpSeq[x].isize() > minDumpLen) {
                    bGood = true;
                    break;
                }
            }
            
            if (!bGood) // no inchworm contig > minDumpLen
                continue;
            
            vector<string> iworm_bundle;
            for (int x = 0; x < tmpSeq.isize(); x++) {
                iworm_bundle.push_back(get_seq_string(tmpSeq[x]));    
                                    
            }
            
            bundled.push_back(iworm_bundle);
            component_ids.push_back(component_no);

            
            //TranscriptomeGraph(tmpSeq,
            //                 pOut,
            //                 24);
            
            tmpSeq.resize(0);
            
        }
        
    }
    
    if (tmpSeq.isize() > 0) {
        
        // capture last one

        if (tmpSeq.isize() > 1 || tmpSeq[0].isize() >= minLen) {
            

            vector<string> iworm_bundle;
            for (int x = 0; x<tmpSeq.isize(); x++) {
                iworm_bundle.push_back(get_seq_string(tmpSeq[x]));
            }


            bundled.push_back(iworm_bundle);
            component_ids.push_back(component_no);
            
                        
            //TranscriptomeGraph(tmpSeq,
            //pOut,
            //                 24);
                        
        }
    }
    
    
    write_iworm_bundle(bundledName, bundled, component_ids, mgr);
    
        
    /////////////////////////////////////////////////////////////////////
    // Make fasta files for each component...
    
    cerr << "-- mapping reads to chrysalis components." << endl;
    
    string readcounts_file = outDir + "/readcounts.out";
    string readsToTranscriptsCompleteFile = outDir + "/readsToTranscripts.finished";
    
    if (! (Exists(readcounts_file) && Exists(readsToTranscriptsCompleteFile) ) ) {

        stringstream command;

        
        
        command << exec << "ReadsToTranscripts -i "
                << readString
                << " -f  "
                << bundledName
                << " -o "
                << outDir;
        if (sStrand) {
            command << " -strand ";
        }
        if (min_pct_read_mapping > 0) {
            command << " -p " << min_pct_read_mapping;
        }

        
        
        if (max_mem_reads > 0) {
            //char max_mem_reads_string[256];
            //sprintf(max_mem_reads_string, " -max_mem_reads %li ",
			//		  max_mem_reads);
			//  command += max_mem_reads_string;
            
            command << " -max_mem_reads " << max_mem_reads;

        }
        //if (bStrand)
        //command += " -strand";
        cout << "Running: " << command.str() << endl;
        Execute(command.str().c_str());
        

        /*

        // Run BWA to map reads to inchworm transcripts
        stringstream cmdstr;
        cmdstr << exec 
               << "/../util/alignReads.pl "
               << " --target " << iwormString
            //<< " --aligner bwa "
               << " --aligner bowtie "
               << " --seqType fa "
               << " --single " << readString
            //<< " -o iworm_bwa "
               << " -o iworm_bowtie --retain_intermediate_files ";
        
        if (sStrand) {
            cmdstr << " --SS_lib_type F ";
        }
        
        // add custom options to aligner
        cmdstr << " -- -a -m 20 --best --strata "; 
        
        Execute(cmdstr.str().c_str());
        
        // distribute the reads according to component
        stringstream cmdstr2;
        if (sStrand) {

            cmdstr2 << exec
                //<< "/../util/distribute_sam_aligned_reads.pl iworm_bwa/iworm_bwa.coordSorted.sam chrysalis/ ";
                    << "/../util/distribute_sam_aligned_reads.pl iworm_bowtie/iworm_bowtie.coordSorted.sam.+.sam chrysalis/ ";
        }
        else {

            cmdstr2 << exec
                //<< "/../util/distribute_sam_aligned_reads.pl iworm_bwa/iworm_bwa.coordSorted.sam chrysalis/ ";
                    << "/../util/distribute_sam_aligned_reads.pl iworm_bowtie/iworm_bowtie.coordSorted.sam chrysalis/ ";
        }
        
        Execute(cmdstr2.str().c_str());
        
        */
        
        
        string complete_cmd = "touch " + readsToTranscriptsCompleteFile;
        Execute(complete_cmd.c_str());
    }
    else {
        cerr << "File: " << readsToTranscriptsCompleteFile << " already exists. Using existing read/transcript mappings assuming resume mode." << endl << endl;
    }
    
    
    // readcounts_file created by the above.
    FlatFileParser readCount;  
    readCount.Open(readcounts_file);
    readCount.ParseLine();
    string numReads = readCount.Line();
        
    
    
    //////////////////////////////////////////////////////////////////////
    // QuantifyGraph
    //////////////////////////////////////////////////////////////////////

    cerr << "-- writing quantifyGraph and butterfly commands for parallel processing." << endl;
    
    //svec<string> targets;
    string butterflyCommands = outDir + "/butterfly_commands";
    FILE * pButterfly = fopen(butterflyCommands.c_str(), "w");

    string quantifyGraphCommands = outDir + "/quantifyGraph_commands";
    FILE * pQGraphCmds = fopen(quantifyGraphCommands.c_str(), "w");

    ofstream component_listing_ofh;
    string component_listing_filename = outDir + "/component_file_listing.txt";
    component_listing_ofh.open(component_listing_filename.c_str());
    
    map<int,bool> pursue_component;

    for (int i=0; i<component_ids.size(); i++) {
        int component_no = component_ids[i];
        string graph = mgr.GetFileName(component_no, ".raw.graph");
        string fasta = mgr.GetFileName(component_no, ".raw.fasta");
        string finalgraph = mgr.GetFileName(component_no, ".out");
        //string fnialreads = mgr.GetFileName(k, ".finalgraph.reads");
        
        string component_file_basename = mgr.GetFileName(component_no, "");
        
        FILE * pFasta = fopen(fasta.c_str(), "r");
        if (pFasta == NULL)
            continue;
       
        fclose(pFasta);
        
        pursue_component[component_no] = true;
        
        char cwd[FILENAME_MAX];
        getcwd(cwd, sizeof(cwd) / sizeof(char));

        stringstream bfly_cmd;

        bfly_cmd << "java -jar " << butterflyExec << " -N " << numReads << " -L " << minLen << " -F " << pairDist;
        bfly_cmd << " -C " << cwd << "/" << mgr.GetFileName(component_no, "");
        
        if (NO_CLEANUP) {
            bfly_cmd << " --no_cleanup ";
        }


        fprintf(pButterfly, "%s\n", bfly_cmd.str().c_str());
        

        stringstream qgraph_cmd;
        qgraph_cmd << exec << "QuantifyGraph -g " << cwd << "/" << graph 
                   << " -o " << cwd << "/" << finalgraph << " -i " << cwd << "/" << fasta;
        if (sStrand) {
            qgraph_cmd << " -strand";
        }
        

        if (max_reads > 0) {
            char max_reads_string[256];
            sprintf(max_reads_string, " -max_reads %li ", max_reads);
            // command += max_reads_string;
            qgraph_cmd << " -max_reads " << max_reads;
        }
        
        if (NO_CLEANUP) {
            qgraph_cmd << " -no_cleanup ";
        }
        
        
        fprintf(pQGraphCmds, "%s\n", qgraph_cmd.str().c_str());

        component_listing_ofh << component_file_basename << endl;




    }
    fclose(pButterfly);
    fclose(pQGraphCmds);
    component_listing_ofh.close();

    // write the deBruijn graphs:
    write_deBruijn_graphs(bundled, component_ids, mgr, outDir, exec, sStrand, pursue_component);
    

    
    return(0);

}
