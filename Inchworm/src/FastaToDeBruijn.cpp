#include <stdlib.h>
#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iterator>
#include <math.h>

#include "Fasta_reader.hpp"
#include "sequenceUtil.hpp"
#include "argProcessor.hpp"
#include "DeBruijnGraph.hpp"
#include "string_util.hpp"
#include "irke_common.hpp"


using namespace std;

unsigned int IRKE_COMMON::MONITOR = 0;





// prototypes
string usage (ArgProcessor args);
int fastaToDeBruijn (int argc, char* argv[]);
DeBruijnGraph constructDeBruijnGraph (vector<string> fasta_filenames, int kmer_length, int component_val, bool sStrand);





string usage (ArgProcessor args) {
    
    stringstream usage_info;
    
    usage_info 
        << endl << endl
        
        << "**Required" << endl
        << "  --fasta  <str>      " << ":fasta file containing reads" << endl
        << "  -K  <int>           " << ":fasta file containing kmers" << endl
        << "  -C <int>            " << ":component identifier" << endl
        << endl
        << " **Optional" << endl
        << "  --SS                " << ":indicates strand-specific" << endl 
        << "  --toString          " << ": dump graph as descriptive output" << endl
        << "  --monitor <int>     " << ": verbosity level" << endl
        << endl;
    
    return(usage_info.str());
    
}



int fastaToDeBruijn (int argc, char* argv[]) {

    try {
        
        ArgProcessor args(argc, argv);
        
        /* Check for essential options */
        if (args.isArgSet("--help") || args.isArgSet("-h")
            ||  (! (args.isArgSet("--fasta") && args.isArgSet("-K") && args.isArgSet("-C")) ) 
            ) {
            
            cerr  << usage(args) << endl << endl;
            
            return 1;
        }
    
        // required
        string fasta_filename = args.getStringVal("--fasta");
        int kmer_length = args.getIntVal("-K");
        int component_val = args.getIntVal("-C");
        
        // optional
        bool sStrand = false;
        if (args.isArgSet("--SS")) {
            sStrand = true;
        }

        if (args.isArgSet("--monitor")) {
            IRKE_COMMON::MONITOR = args.getIntVal("--monitor");
        }
        
        vector<string> fasta_filenames;
        if (fasta_filename.find(',') != string::npos) {
            // cerr << "Splitting filenames out: " << fasta_filename << endl;
            string_util::tokenize(fasta_filename, fasta_filenames, ",");
        }
        else {
            //cerr << "Using single fasta filename: " << fasta_filename << endl;
            fasta_filenames.push_back(fasta_filename);
        }
        
        
        DeBruijnGraph g = constructDeBruijnGraph(fasta_filenames, kmer_length, component_val, sStrand);
        
        if (args.isArgSet("--toString")) {
            cout << g.toString();
        }
        else {
            cout << g.toChrysalisFormat(component_val, sStrand);
        }
        
    }
    
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    
    
    
    return(0);
    
}

DeBruijnGraph constructDeBruijnGraph (vector<string> fasta_file_names, int kmer_length, int component_val, bool sStrand) {


    DeBruijnGraph g(kmer_length);
    
    for (int i = 0; i < fasta_file_names.size(); i++) {
        
        string fasta_filename = fasta_file_names[i];
        
        if (IRKE_COMMON::MONITOR > 1)
            cerr << "Parsing file: " << fasta_filename << endl;
        
        Fasta_reader fasta_reader(fasta_filename);
        
        
        
        while (fasta_reader.hasNext()) {
            
            Fasta_entry fe = fasta_reader.getNext();
            string sequence = fe.get_sequence();
            
            vector<string> seq_regions;
            string_util::tokenize(sequence, seq_regions, "X"); // inchworm bundles concatenated with 'X' delimiters by Chrysalis 
            
            for (int s = 0; s < seq_regions.size(); s++) { 
                
                string seq_region = seq_regions[s];
                
                if (contains_non_gatc(seq_region)) {
                    seq_region = replace_nonGATC_chars_with_A(seq_region);
                }
                

                if (IRKE_COMMON::MONITOR > 2) 
                    cerr << "Adding sequence to graph: " << seq_region << endl;

                g.add_sequence(seq_region);
                
                if (! sStrand) {
                    string revseq = revcomp(seq_region);
                    
                    if (IRKE_COMMON::MONITOR > 2) 
                        cerr << "Adding sequence to graph: " << revseq << endl;


                    g.add_sequence(revseq);
                    
                    

                }
            }
        }
    }
    
    return(g);
    
    
    

} 



    
int main (int argc, char* argv[]) {
    
    try {
        return(fastaToDeBruijn(argc, argv));
    }
    
    catch (string err) {
        cerr << err << endl;
    }
    
    return(1);
    
}


