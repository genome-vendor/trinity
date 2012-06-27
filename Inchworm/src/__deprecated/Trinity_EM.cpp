#include <stdlib.h>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <iterator>
#include <iostream>
#include <math.h>
#include <fstream>
#include "argProcessor.hpp"
#include "SAM_entry.hpp"
#include "SAM_reader.hpp"
#include "Fasta_reader.hpp"
#include "Fasta_entry.hpp"
#include "string_util.hpp"
#include "EM.hpp"
#include "common.hpp"

using namespace std;


unsigned int DEBUG_LEVEL = 0;


// function protos
int run_Trinity_EM (int argc, char* argv[]);

string usage(ArgProcessor args);

void process_cmd_line_arguments(int argc, char* argv[]);


string extract_component_from_transcript_name (const string& transcript);


map<string,bool> examine_EM_per_component_phase_1(const string& coord_sorted_sam_file, const map<string,unsigned int>& transcript_seq_len_map, const unsigned int fragment_length);

void run_component_EM_phase_1 (map<string, map<string,bool> >& read_to_trans, 
                               map<string,bool>& transcripts_OK, 
                               const map<string,unsigned int>& transcript_seq_len_map, 
                               const string& curr_component,
                               const unsigned int fragment_length);

map<string,bool> identify_repetitive_reads(const string& coord_sorted_sam_file);

map<string, map<string,bool> >& require_unique_read_mapping (map<string, map<string,bool> >& read_to_trans, unsigned int min_unique_read_mappings_required);

double compute_sum_fpkm(const vector<t_EM_result>& em_results);

void compute_abundance_estimates_filtered_Trinity_phase_2(map<string,bool>& transcripts_passed, 
                                                          const string& NAME_SORTED_SAM_FILE, 
                                                          const map<string,unsigned int>& transcript_seq_len_map,
                                                          const map<string,bool>& repetitive_reads,
                                                          const unsigned int fragment_length);


bool examine_components_for_filtering_phase_2(const vector<t_EM_result>& em_results, map<string,bool>& transcripts_OK);




map<string, vector<t_EM_result> > group_EM_results_by_component(const vector<t_EM_result>& em_results);

map<string,bool> get_transcript_map_from_seq_map (const map<string,unsigned int>& transcript_seq_len_map);

map<string,unsigned int> get_trans_lengths(const map<string,string>&  trans_seq_map);



/* vars to store command-line parameters and defaults */
// required arguments.
string COORD_SORTED_SAM_FILE = "";
string NAME_SORTED_SAM_FILE = "";
string TRINITY_FASTA_FILE = "";
unsigned int FRAGMENT_LENGTH = 0;
string OUT_PREFIX = "Trinity.refined";

// phase 1
float RETAIN_MIN_PERCENT_COMPONENT_EXPRESSED = 1.0f;
int MAX_ITERATIONS_PHASE_1 = 1000;
int ITERATIONS_REVIEW_STEP_PHASE_1 = 1;
bool AGGRESSIVE = false;
unsigned int MAX_REPETITIVE_READ = 3;
int MIN_UNIQUE_READ_PHASE_1 = 0;


// phase 2
int MAX_ITERATIONS_PHASE_2 = 1000;

// both phases
float MIN_DELTA_ML = 0.01f;


bool ONE_COMPONENT_MODE = false;


// SAM parsing settings
bool REQUIRE_PROPERLY_MAPPED_PAIRS = 0;

int main (int argc, char* argv[]) {

    try {
        return(run_Trinity_EM(argc, argv));
    }
    catch (string err) {
        cerr << err << endl;
    }

    return(1);
}


int run_Trinity_EM (int argc, char* argv[]) {

    process_cmd_line_arguments(argc, argv);

    Fasta_reader fasta_reader(TRINITY_FASTA_FILE);
    map<string,string> transcript_seqs = fasta_reader.retrieve_all_seqs_hash();
    map<string,unsigned int> transcript_seq_len_map = get_trans_lengths(transcript_seqs);
    
    // examine each component's data ata time
   
    map<string,bool> transcripts_passed;
    
    
    if (! ONE_COMPONENT_MODE) {
        // Phase 1: examine expression on per-component basis.
        cerr << "\n\n\n#######   Phase 1, Per-Component Analysis/Filtering ##########\n\n";
        transcripts_passed = examine_EM_per_component_phase_1(COORD_SORTED_SAM_FILE, transcript_seq_len_map, FRAGMENT_LENGTH);
    }
    else {
        transcripts_passed = get_transcript_map_from_seq_map(transcript_seq_len_map);

    }

    map<string,bool> repetitive_reads;
    
    if (MAX_REPETITIVE_READ > 1) {
        cerr << "\n\n## Identifying repetitive reads" << endl;
        repetitive_reads = identify_repetitive_reads(COORD_SORTED_SAM_FILE);
    }

    
    // Phase 2: report final transcripts and abundance estimates.
    cerr << "\n\n\n########   Phase 2, Analysis Exploring All Reads and All Components Simutaneously  ###########\n\n\n";

    compute_abundance_estimates_filtered_Trinity_phase_2(transcripts_passed, NAME_SORTED_SAM_FILE, transcript_seq_len_map, repetitive_reads, FRAGMENT_LENGTH);
    
    cerr << "\n\nFinished.\n\n";
    
    
    return(0);
}
    


void process_cmd_line_arguments (int argc, char* argv[]) {

    try {
        ArgProcessor args(argc, argv);

        /* check for essential options */
        if (args.isArgSet("--help") || args.isArgSet("-h")
            ||
            (! 
             (args.isArgSet("--name_sorted_sam")
              &&
              args.isArgSet("--fasta") 
              &&
              args.isArgSet("--fragment_length")
              
              )
             )
            )
            {
                cerr << usage(args) << endl << endl;
                exit(1);
            }
        
        if (args.isArgSet("--debug_level")) {
            DEBUG_LEVEL = args.getIntVal("--debug_level");
        }
        
        if (args.isArgSet("--coord_sorted_sam")) {
            COORD_SORTED_SAM_FILE = args.getStringVal("--coord_sorted_sam");
        }
        if (args.isArgSet("--name_sorted_sam")) {
            NAME_SORTED_SAM_FILE = args.getStringVal("--name_sorted_sam");
        }
        if (args.isArgSet("--fasta")) {
            TRINITY_FASTA_FILE = args.getStringVal("--fasta");
        }

        if (args.isArgSet("--out_prefix")) {
            OUT_PREFIX = args.getStringVal("--out_prefix");
        }
        
        if (args.isArgSet("--retain_min_percent_component_expressed")) {
            RETAIN_MIN_PERCENT_COMPONENT_EXPRESSED = args.getFloatVal("--retain_min_percent_component_expressed");
        }
        if (args.isArgSet("--max_iterations_phase_1")) {
            MAX_ITERATIONS_PHASE_1 = args.getIntVal("--max_iterations_phase_1");
        }
        if (args.isArgSet("--min_unique_read_phase_1")) {
            MIN_UNIQUE_READ_PHASE_1 = args.getIntVal("--min_unique_read_phase_1");
        }
                
        if (args.isArgSet("--aggressive")) {
            AGGRESSIVE = true;
        }
        if (args.isArgSet("--max_iterations_phase_2")) {
            MAX_ITERATIONS_PHASE_2 = args.getIntVal("--max_iterations_phase_2");
        }
        if (args.isArgSet("--min_delta_ML")) {
            MIN_DELTA_ML = args.getFloatVal("--min_delta_ML");
        }
        if (args.isArgSet("--max_repetitive_read")) {
            MAX_REPETITIVE_READ = args.getIntVal("--max_repetitive_read");
        }

        if (args.isArgSet("--one_component_mode")) {
            ONE_COMPONENT_MODE = true;
        }


        if (args.isArgSet("--require_properly_mapped_pairs")) {
            REQUIRE_PROPERLY_MAPPED_PAIRS = true;
        }
        
        if (args.isArgSet("--fragment_length")) {
            FRAGMENT_LENGTH = args.getIntVal("--fragment_length");
        }
        
        
        if (MAX_REPETITIVE_READ > 0 || ! ONE_COMPONENT_MODE) {
            if (! args.isArgSet("--coord_sorted_sam")) {
                cerr << "If --max_repetitive_read > 0, you must set --coord_sorted_sam" << endl;
                exit(2);
            }
        }
        

    }
    
    catch (exception& e) {
        cerr << "error: " << e.what() << "\n";
        exit(1);
    }
    

    return;
    
   
}


string usage (ArgProcessor args) {
    
    stringstream usage_info;

    usage_info
        << endl << endl
        << "#################################################################################################################" << endl
        << "#" << endl
        << "#  ## Required:" << endl
        << "#" << endl
        << "#      --fasta <string>                   Trinity.fasta" << endl        
        << "#      --coord_sorted_sam <string>        coord-sorted sam  (  sort -k3,3 -k4,4n file.sam > coordSorted.sam )" << endl
        << "#      --name_sorted_sam <string>         read name-sorted sam  (   sort -k1,1 -k3,3 file.sam > nameSorted.sam )" << endl
        << "#      --fragment_length <int>            mean length of a fragment from which a RNA-Seq read was derived." << endl
    
        << "#" << endl
        << "#  ## Optional" << endl
        << "#" << endl
        << "#      --out_prefix <string>              prefix for output files. (default: " << OUT_PREFIX << ", resulting in " << endl
        << "#                                              output file: " << OUT_PREFIX << "..expression)" << endl
        << "#" << endl
        << "#      --retain_min_percent_component_expressed <float>  (default: " << RETAIN_MIN_PERCENT_COMPONENT_EXPRESSED << ")" << endl
        << "#" << endl
        << "#  * Phase 1 (component pre-screen):" << endl
        << "#      --max_iterations_phase_1  <int>    component-based pre-screen, maximal number of EM iterations (default: " << MAX_ITERATIONS_PHASE_1 << ")" << endl
        // dangerous, potential to remove transcripts that are best
        // needs re-thinking here
        // << "#      --min_unique_read_phase_1 <int>    minimum number of unique reads required to belong to a transcript within a " << endl
        // << "#                                              component (default: " << MIN_UNIQUE_READ_PHASE_1 << ") " << endl
        // << "#      --aggressive                       rapidly removes transcripts that have no unique read mapping.  By default, " << endl
        // << "#                                              deletes one at a time and reexamines uniqueness after each deletion." << endl
        << "#      --max_repetitive_read <int>         maximum number of components that a single read can map to before being " << endl
        << "#                                              eliminated (default: " << MAX_REPETITIVE_READ << ")" << endl
        << "#      --one_component_mode                treat entire data set as one large component rather than the individual " << endl
        << "#                                              Trinity components (skips phase 1 altogether, and can be extremely slow)" << endl
        << "#" << endl
        << "#  * Phase 2: (full data set simultaneously)" << endl
        << "#      --max_iterations_phase_2 <int>     final estimation of abundance, max number of EM interations (default: " << MAX_ITERATIONS_PHASE_2 << ")" << endl
        << "#" << endl
        << "#  Both phases" << endl 
        << "#      --min_delta_ML                      EM stops when difference between subsequent max likelihood computations " << endl
        << "#                                              is below this value. (default: " << MIN_DELTA_ML << ")" << endl
        << "#" << endl
        << "#  * SAM parsing" << endl
        << "#      --require_properly_mapped_pairs     In the context of paired reads, only considers mappings involving properly-" << endl
        << "#" << endl
        << "   * misc" << endl
        << "#      --debug_level <int> " << endl
        << "#                                              mapped pairs. " << endl
        << "#" << endl
        << "#" << endl
        << "#####################################################################################################################" << endl << endl;

    
    return(usage_info.str());
}


map<string,bool> examine_EM_per_component_phase_1(const string& coord_sorted_sam_file, const map<string,unsigned int>& transcript_seq_len_map, const unsigned int fragment_length) {
    
    map<string,bool> transcripts_OK;

    string curr_component = "";

    map<string, map<string,bool> > read_to_trans;

    cerr << "-processing coord_sorted_sam_file: " << coord_sorted_sam_file << endl;
    
    unsigned long read_counter = 0;
    
    SAM_reader sam_reader(coord_sorted_sam_file);

    while (sam_reader.has_next()) {
        read_counter++;
        SAM_entry sam_entry = sam_reader.get_next();
        
        if (REQUIRE_PROPERLY_MAPPED_PAIRS) {
            if (sam_entry.is_paired() && ! sam_entry.is_proper_pair()) {
                continue;
            }
        }
        

        if (read_counter % 100000 == 0) {
            cerr << "\r[" << read_counter << "] sam entries read     ";
        }
        
        string trans_id = sam_entry.get_scaffold_name();
        string read_name = sam_entry.get_read_name();
        string component_id = extract_component_from_transcript_name(trans_id);
        
        if (component_id != curr_component) {
            if (! read_to_trans.empty()) {
                run_component_EM_phase_1(read_to_trans, transcripts_OK, transcript_seq_len_map, curr_component, fragment_length);
                read_to_trans.clear(); // reinit
            }
        }
        curr_component = component_id;
        if (read_to_trans.find(read_name) != read_to_trans.end()) {
            map<string, map<string,bool> >::iterator it  = read_to_trans.find(read_name);
            it->second[trans_id] = true;
        }
        else {
            // insert it
            map<string,bool> read_map;
            read_map[trans_id] = true;
            read_to_trans[read_name] = read_map;
            
        }
        
        
    }
    
    if (! read_to_trans.empty()) {
        run_component_EM_phase_1(read_to_trans, transcripts_OK, transcript_seq_len_map, curr_component, fragment_length);
    }
    
    return(transcripts_OK);
}


void run_component_EM_phase_1 (map<string, map<string,bool> >& read_to_trans,
                               map<string,bool>& transcripts_OK,
                               const map<string,unsigned int>& transcript_seq_len_map,
                               const string& curr_component,
                               const unsigned int fragment_length) {
    
    
    EM em_obj(transcript_seq_len_map, fragment_length);

    if (MIN_UNIQUE_READ_PHASE_1 > 0) {
        cerr << endl << "-removing transcripts having less than " << MIN_UNIQUE_READ_PHASE_1 << endl;
        read_to_trans = require_unique_read_mapping(read_to_trans, MIN_UNIQUE_READ_PHASE_1);
    }
    
    // add reads
    for (map<string, map<string,bool> >::iterator it = read_to_trans.begin();
         it != read_to_trans.end();
         it++) {

        string read_name = it->first;
        map<string,bool> transcripts_mapped = it->second;
        
        vector<string> transcripts;
        for (map<string,bool>::iterator it2 = transcripts_mapped.begin(); it2 != transcripts_mapped.end(); it2++) {
            
            string transcript = it2->first;
            transcripts.push_back(transcript);
        }
        
        em_obj.add_read(transcripts);
        
	}
        
    
    
    double prev_ML = em_obj.run(10, MIN_DELTA_ML);;
    
    for (int i = 11; i <= MAX_ITERATIONS_PHASE_1; i++) {
        
        double ML = em_obj.resume(1, MIN_DELTA_ML, prev_ML);
        em_obj.purge_low_percent_expressed_transcripts(RETAIN_MIN_PERCENT_COMPONENT_EXPRESSED);
        
        if (ML - prev_ML <= MIN_DELTA_ML) {
            break;
        }
        prev_ML = ML;
    }
    cerr << endl;
    
    
    vector<t_EM_result> em_results = em_obj.get_and_report_results();
    

    // filter those transcripts that are lowly expressed
    double sum_fpkm = compute_sum_fpkm(em_results);

    for (vector<t_EM_result>::iterator it = em_results.begin(); it != em_results.end(); it++) {
        
        t_EM_result result = *it;
        
        double fpkm = result.FPKM;

        double percent_component_expressed = fpkm/sum_fpkm * 100;

        if (percent_component_expressed >= RETAIN_MIN_PERCENT_COMPONENT_EXPRESSED) {
            
            string transcript = result.trans_id;
            
            transcripts_OK[ transcript ] = true;
        }
    
    }
    
    return;

    
    
}




string extract_component_from_transcript_name (const string& transcript) {
    
    if (ONE_COMPONENT_MODE) {
        return("ONE_COMPONENT");
    }
    else {

        vector<string> tokens;
        string_util::tokenize(transcript, tokens, "_");
        
        return(tokens[0]);
    }
}





map<string, map<string,bool> >& require_unique_read_mapping (map<string, map<string,bool> >& read_to_trans, unsigned int min_unique_read_mappings_required) {
    
    // only do aggressive mapping for now.
    // TODO: update for less aggressive deletion.
    
    
    // capture uniquely mapped transcripts.
    
    int round = 0;

    while (true) {

        round++;
        cerr << "-round[" << round << "] removing low support trans. " << endl;
        
        map<string, unsigned long> transcript_to_read_abundance;
        map<string, unsigned long> transcript_to_unique_read_count;
        
        for (map<string, map<string,bool> >::iterator it = read_to_trans.begin();
             it != read_to_trans.end();
             it++) {
            
            string read_name = it->first;
            
            map<string,bool> transcript_mapping = it->second;
            
            if (transcript_mapping.size() == 1) {
                
                string transcript = transcript_mapping.begin()->first;
                transcript_to_unique_read_count[ transcript ] ++;
                transcript_to_read_abundance[ transcript ] ++;
            }
            else {
                for (map<string,bool>::iterator it2 = transcript_mapping.begin(); it2 != transcript_mapping.end(); it2++) {
                    string transcript = it2->first;
                    transcript_to_read_abundance[ transcript ] ++;
                }
            }
        }
        
        
        // Summarize all data for Round 1
        if (round == 1) {
            for (map<string, unsigned long>::iterator it = transcript_to_read_abundance.begin();
                 it != transcript_to_read_abundance.end();
                 it++) { 
                
                string transcript = it->first;
                unsigned long read_abundance = it->second;
                
                unsigned long unique_read_count = transcript_to_unique_read_count[transcript];
                
                cerr << "Initial: " << transcript << " unique: " << unique_read_count << " abundance: " << read_abundance << endl;
            }
        }
        

        // capture transcripts that lack unique reads.
        
        map<string,bool> insufficient_unique_read_mapped_transcripts;
        for (map<string, unsigned long>::iterator it = transcript_to_read_abundance.begin(); it != transcript_to_read_abundance.end(); it++) {
            
            string transcript = it->first;
            
            map<string, unsigned long>::iterator it2 = transcript_to_unique_read_count.find(transcript);
            if (it2 == transcript_to_unique_read_count.end() || it2->second < min_unique_read_mappings_required) {
                insufficient_unique_read_mapped_transcripts[transcript] = true;
            }
            
        }
        
        if (insufficient_unique_read_mapped_transcripts.empty()) {
            break; // done filtering
        }
        
                
        if (AGGRESSIVE) {
            // remove all such transcripts.
            cerr << "-AGGRESSIVE mode, removing " << insufficient_unique_read_mapped_transcripts.size() << " transcripts with insufficient support." << endl;
        }
        else {
            // remove the one transcript with the smallest number of uniquely mapped reads followed by the smallest number of multiply mapped reads
            string transcript_to_remove;
            unsigned long min_unique_read_count = min_unique_read_mappings_required;
            unsigned long min_abundance = 100000; // just a high number, will get replaced shortly
            for (map<string, bool>::iterator it = insufficient_unique_read_mapped_transcripts.begin();
                 it != insufficient_unique_read_mapped_transcripts.end();
                 it++) {
                
                string transcript = it->first;
                unsigned long uniq_count = transcript_to_unique_read_count[transcript];
                unsigned long abundance = transcript_to_read_abundance[transcript];
                
                if (uniq_count < min_unique_read_count ||
                    (uniq_count == min_unique_read_count && min_abundance > abundance)
                    ) {
                    transcript_to_remove = transcript;
                    min_unique_read_count = uniq_count;
                    min_abundance = abundance;
                }
            }
            
            // replace with the one transcript to purge
            insufficient_unique_read_mapped_transcripts.clear();
            insufficient_unique_read_mapped_transcripts[transcript_to_remove] = true;
            cerr << "   (removing: " << transcript_to_remove << ", uniq_count: " << min_unique_read_count << ", total read abundance: " << min_abundance << ")" << endl;
        }
        
        // iterate read -> transcripts, remove offending transcripts from read mappings
        vector<string> reads_to_purge; // those reads that now have no mappings
        for (map<string, map<string,bool> >::iterator it = read_to_trans.begin();
             it != read_to_trans.end();
             it++) {
            
            string read_name = it->first;
            
            map<string,bool>& transcript_mapping = it->second;
            
            map<string,bool> scheduled_deletion;
            
            for (map<string,bool>::iterator it2 = transcript_mapping.begin(); it2 != transcript_mapping.end(); it2++) {
                string transcript = it2->first;
                if (insufficient_unique_read_mapped_transcripts.find(transcript) != insufficient_unique_read_mapped_transcripts.end()) {
                    scheduled_deletion[transcript] = 1;
                }
                
            }
            
            for (map<string,bool>::iterator it3 = scheduled_deletion.begin(); it3 != scheduled_deletion.end(); it3++) {
                string transcript = it3->first;
                transcript_mapping.erase(transcript);
            }
            
            if (transcript_mapping.empty()) {
                reads_to_purge.push_back(read_name);
            }
            
                
            
        }

        for (vector<string>::iterator it = reads_to_purge.begin(); it != reads_to_purge.end(); it++) {
            string read_to_delete = *it;
            read_to_trans.erase(read_to_delete);
        }
        

    }
    
    return(read_to_trans);
}


double compute_sum_fpkm(const vector<t_EM_result>& em_results) {

    double sum_fpkm = 0.0;

    for (vector<t_EM_result>::const_iterator it = em_results.begin(); it != em_results.end(); it++) {
        
        t_EM_result result = *it;
        
        sum_fpkm += result.FPKM;
    }

    return(sum_fpkm);
}


void compute_abundance_estimates_filtered_Trinity_phase_2(map<string,bool>& transcripts_OK,
                                                          const string& NAME_SORTED_SAM_FILE,
                                                          const map<string,unsigned int>& transcript_seq_len_map,
                                                          const map<string,bool>& repetitive_reads,
                                                          const unsigned int fragment_length) {
    
    
    cerr << "-there are " << transcripts_OK.size() << " transcripts that pass minimal expression criteria." << endl;
    
    cerr << "-processing name_sorted_sam_file: " << NAME_SORTED_SAM_FILE << endl;

    bool filtered_transcripts_flag = 1;
    
    vector<t_EM_result> em_results;  // update each round.
    unsigned long total_mapped_read_count;
    unsigned long total_num_transcripts;

    int round = 0;

    while (filtered_transcripts_flag) {
        round++;
        
        cerr << "-phase 2, round: " << round << endl;

        EM em_obj (transcript_seq_len_map, fragment_length); // reinit each round

        SAM_reader sam_reader(NAME_SORTED_SAM_FILE);
        
        unsigned long read_counter = 0;
        unsigned long repetitive_read_counter = 0;
        unsigned long sam_line_counter = 0;
        
        string curr_read_name = "";
        map<string,bool> transcripts;

        while (sam_reader.has_next()) {
            
            SAM_entry sam_entry = sam_reader.get_next();

            sam_line_counter++;
            
            if (REQUIRE_PROPERLY_MAPPED_PAIRS) {
                if (sam_entry.is_paired() && ! sam_entry.is_proper_pair()) {
                    continue;
                }
            }

            string read_name = sam_entry.get_read_name();
            string transcript = sam_entry.get_scaffold_name();

            if (repetitive_reads.find(read_name) != repetitive_reads.end()) {
                // got a repetitive read
                repetitive_read_counter++;
                continue;
            }
            

            if (transcripts_OK.find(transcript) == transcripts_OK.end()) {
                // not OK
                continue;
            }

            
            
            if (curr_read_name != "" && curr_read_name != read_name) {

                // process current collection of transcripts corresponding to the previously tracked read.
                
                read_counter++;
                if (read_counter % 100000 == 0) {
                    cerr << "\r[" << read_counter << "] reads read, " 
                         << sam_line_counter << " sam lines, of which " 
                         << repetitive_read_counter << " were ignored as repetitive.";
                }
                if (! transcripts.empty()) {
                    vector<string> transcript_list;
                    for (map<string,bool>::iterator it = transcripts.begin(); it != transcripts.end(); it++) {
                        transcript_list.push_back(it->first);
                    }
                    em_obj.add_read(transcript_list);
                }
                transcripts.clear();
            }

            transcripts[transcript] = true;
                        
            curr_read_name = read_name;
        }

        // get last ones
        if (! transcripts.empty()) {
            vector<string> transcript_list;
            for (map<string,bool>::iterator it = transcripts.begin(); it != transcripts.end(); it++) {
                transcript_list.push_back(it->first);
            }
            em_obj.add_read(transcript_list);
            
        }
        
        total_mapped_read_count = em_obj.get_total_mapped_read_count();
        total_num_transcripts = em_obj.get_num_transcripts();
        
        cerr << "-done parsing reads, now running EM" << endl;
        
        em_obj.run(MAX_ITERATIONS_PHASE_2, MIN_DELTA_ML);
        
        em_results = em_obj.get_and_report_results();
        
        if (ONE_COMPONENT_MODE) {
            // no filtering
            break;
        }
        
        filtered_transcripts_flag = examine_components_for_filtering_phase_2(em_results, transcripts_OK);
        
    }
    
    cerr << "Process is now complete. Now writing result files...." << endl;

    // Report final results:
    
    stringstream expression_report;
    
    expression_report << "#Total reads mapped: " << total_mapped_read_count << endl;
    
    expression_report << "#Total transcripts examined: " << total_num_transcripts<< endl;

    expression_report << "#transcript" << "\ttrans_length" <<  "\tunique_map" << "\tmulti_map"
                      << "\tEM_frag_count" << "\tFPKM" <<  "\t%ComponentExpressed" << endl;
    
    map<string, vector<t_EM_result> > component_to_trans_results = group_EM_results_by_component(em_results);


    // write output files
    string expression_report_filename = OUT_PREFIX + ".expression_report";
    ofstream expr_ofh (expression_report_filename.c_str());
        
    for (map<string, vector<t_EM_result> >::iterator it = component_to_trans_results.begin();
         it != component_to_trans_results.end();
         it++) {
        
        string component_id  = it->first;
        
        vector<t_EM_result> component_em_results = it->second;
        
        double sum_fpkm = compute_sum_fpkm(component_em_results);
        
        for (vector<t_EM_result>::iterator it2 = component_em_results.begin();
             it2 != component_em_results.end();
             it2++) {
            
            t_EM_result result = *it2;
            
            double percent_component_expressed = result.FPKM/sum_fpkm * 100;
            
            expression_report << result.trans_id << "\t"
                              << result.length << "\t"
                              << result.unique_map << "\t"
                              << result.multi_map << "\t"
                              << result.expected_map << "\t"
                              << result.FPKM << "\t"
                              << percent_component_expressed << "%\n";
            
            
                        
        }
        expression_report << endl; // separator between components.
    }

    cout << expression_report.str();
    expr_ofh << expression_report.str();
    
    expr_ofh.close();
    

    return;
    

}


bool examine_components_for_filtering_phase_2(const vector<t_EM_result>& em_results, map<string,bool>& transcripts_OK) {

    // first group results by component

    bool filtered_transcript_flag = false;
    
    map<string, vector<t_EM_result> > component_to_trans_results = group_EM_results_by_component(em_results);
    
    
    
    // iterate through each component and validate expression requirements
    for (map<string, vector<t_EM_result> >::iterator it = component_to_trans_results.begin();
         it != component_to_trans_results.end();
         it++) {

        string component_id  = it->first;
        
        vector<t_EM_result> component_em_results = it->second;
        
        double sum_fpkm = compute_sum_fpkm(component_em_results);
        
        for (vector<t_EM_result>::iterator it2 = component_em_results.begin();
             it2 != component_em_results.end();
             it2++) {

            t_EM_result result = *it2;
            
            double percent_component_expressed = result.FPKM/sum_fpkm * 100;

            if (percent_component_expressed < RETAIN_MIN_PERCENT_COMPONENT_EXPRESSED) {
                string transcript = result.trans_id;
                transcripts_OK.erase(transcript);
                
                filtered_transcript_flag = true;
                cerr << "-purging transcript: " << transcript << " due to low support: " << percent_component_expressed << "%" << endl;
            }
            
        }
    }
    


    return(filtered_transcript_flag);
}



map<string, vector<t_EM_result> > group_EM_results_by_component(const vector<t_EM_result>& em_results) {

    map<string, vector<t_EM_result> > component_to_trans_results;
    
    for (vector<t_EM_result>::const_iterator it = em_results.begin(); it != em_results.end(); it++) {
        
        t_EM_result result = *it;
        
        string transcript = result.trans_id;
        
        string component_id = extract_component_from_transcript_name(transcript);

        component_to_trans_results[component_id].push_back(result);
                    
        

    }

    return(component_to_trans_results);
    


}


map<string,bool> get_transcript_map_from_seq_map (const map<string,unsigned int>& transcript_seq_len_map) {

    map<string,bool> transcript_accs;

    for (map<string,unsigned int>::const_iterator it = transcript_seq_len_map.begin(); it != transcript_seq_len_map.end(); it++) {
        
        string transcript = it->first;
        transcript_accs[transcript] = true;
    }
    
    return(transcript_accs);
}


map<string,unsigned int> get_trans_lengths(const map<string,string>&  trans_seq_map) {


    map<string,unsigned int> trans_seq_len_map;


    map<string,string>::const_iterator it = trans_seq_map.begin();
    
    
    while (it != trans_seq_map.end()) {
        
        string transcript = it->first;
        string sequence = it->second;
        
        trans_seq_len_map[transcript] = sequence.length();
        
        it++;
    }
     
    return (trans_seq_len_map);

}

map<string,bool> identify_repetitive_reads(const string& coord_sorted_sam_file) {

    map<string, map<string,bool> > read_to_component_map;
    
    SAM_reader sam_reader (coord_sorted_sam_file);
    
    unsigned long sam_read_counter = 0;

    while(sam_reader.has_next()) {
        SAM_entry sam_entry = sam_reader.get_next();
        
        if (REQUIRE_PROPERLY_MAPPED_PAIRS) {
            if (sam_entry.is_paired() && ! sam_entry.is_proper_pair()) {
                continue;
            }
        }


        string read_name = sam_entry.get_read_name();
        string transcript = sam_entry.get_scaffold_name();
        string component_id = extract_component_from_transcript_name(transcript);
        
        read_to_component_map[read_name][component_id] = true;
    
        sam_read_counter++;
        if (sam_read_counter % 1000 == 0) {
            cerr << "\r[" << sam_read_counter << "] sam entries read.     ";
        }
    }
    

    // identify those reads that map to too many components
    map<string,bool> repetitive_reads;
    for (map<string, map<string,bool> >::iterator it = read_to_component_map.begin();
         it != read_to_component_map.end();
         it++) {
        
        string read_name = it->first;
        map<string,bool>& component_map = it->second;
        
        if (component_map.size() >= MAX_REPETITIVE_READ) {
            repetitive_reads[read_name] = true;
        }
    }
    
    
    return(repetitive_reads);
    
}

