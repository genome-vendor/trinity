#include <stdlib.h>
#include "Fasta_reader.hpp"
#include "Fasta_entry.hpp"
#include "string_util.hpp"
#include "sequenceUtil.hpp"
#include "argProcessor.hpp"
#include <string>
#include <iostream>
#include <sstream>
#include "Ktree.hpp"


int run(int argc, char* argv[]);


int main (int argc, char* argv[]) {

    try {
        return(run(argc, argv));
    }
    
    catch (string errmsg) {
        cerr << errmsg;
        return(1);
    }
    
    
}

int run (int argc, char* argv[]) {
    
    if (argc < 3) {
        stringstream s;
        s << "Usage: " << argv[0] << " file.fasta kmer_length [DS_mode]" << endl << endl;
        
        cerr << s.str();
        return(1);
        
    }

    string fasta_filename (argv[1]);
    unsigned int kmer_length = atoi(argv[2]);
    
    bool DS_mode = (argc >= 3) ? true : false;
    
    Fasta_reader fasta_reader(fasta_filename);
    
    Ktree ktree;

    long read_counter = 0;
    
    while (fasta_reader.hasNext()) {
        
        read_counter++;
        if (read_counter % 1000 == 0) {
            cerr << "\rread[" << read_counter << "]   ";
        }
        

        Fasta_entry fe = fasta_reader.getNext();
        
        string accession = fe.get_accession();
        
        
        string sequence = fe.get_sequence();
        
        // cerr << "Processing: " << sequence << endl;
                        
        if (sequence.length() < kmer_length + 1) {
            continue;
        }
        
        for (unsigned int i = 0; i <= sequence.length() - kmer_length; i++) {
            
            string kmer = sequence.substr(i, kmer_length);
            
            if (! contains_non_gatc(kmer)) {

                ktree.add_kmer(kmer);
            
                if (DS_mode) {
                    kmer = revcomp(kmer);
                    ktree.add_kmer(kmer);
                }

            }
            
        }
        
    }
 

    ktree.report_kmer_counts();
    
   
    return(0);
}

