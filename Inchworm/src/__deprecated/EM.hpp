#ifndef __EM__
#define __EM__


#include <map>
#include <string>
#include <vector>
#include <sstream>
#include "stacktrace.hpp"
#include <math.h>
#include "string_util.hpp"

using namespace std;


typedef struct {
    string trans_id;
    unsigned int length;
    unsigned long unique_map;
    unsigned long multi_map;
    double expected_map;
    double theta;
    double FPKM;
} t_EM_result;


typedef struct {
    unsigned long uniquely_mapped_read_count;
    unsigned long multiply_mapped_read_count;
} t_read_counts;


class EM {

public:

    EM(const map<string,unsigned int>& trans_seq_len_map, unsigned int frag_length); 
    
    double run (unsigned int max_iterations, double min_delta_ML); // returns ML value from last round

    double resume (unsigned int max_iterations, double min_delta_ML, double prev_ML_val); // same as above, but run() must have been executed first.
    

    vector<string> get_all_transcripts();
    
    unsigned long get_num_transcripts();
    

    vector<t_EM_result> get_results();
    
    vector<t_EM_result> get_and_report_results();
    
    void add_read (vector<string>& transcripts);    

    t_read_counts count_reads_mapped_to_transcript (const string& transcript);

    
    unsigned long get_total_mapped_read_count();
    
    bool purge_low_percent_expressed_transcripts(float min_percent_expressed);
    

private:

    /*  ---- PRIVATE VARS ----- */
    
    // maping info:   transcript => (transcript_combo => count)  
    //      where count = num reads that map to all in combo, for which transcript is one entry.
    map<string, map<string, unsigned long> > _trans_to_multi_map_counts;
    
    // model params
    map<string, double> _ENt; // expected read count for transcript
    
    map<string, double> _theta; // fraction of all reads corresponding to transcript

    // transcript info
    const map<string, unsigned int>& _trans_lengths;  // lengths of transcripts
    
    const unsigned int _frag_length; // length of a single RNA-Seq fragment from which the read was derived.

    // misc
    unsigned long _total_reads; 
    
    map<string, vector<string> > _combo_memoized;
    
    /* ---------- PRIVATE METHODS  ----------*/
    

    // EM-functions
    void _init_theta();
    void _E_step();
    void _M_step();
    double _compute_Max_Likelihood();
    double _get_theta(const string& transcript);
    double _get_ENt(const string& transcript);
    unsigned int _em_round;
    
    // TODO: rename these methods for improved consistency in use of terms.
    
    map<string, unsigned long>& _get_combo_map_for_transcript (const string& transcript);
    map<string, unsigned long>& _get_or_create_combo_map_for_transcript (const string& transcript);
    
    
    vector<string> _get_multi_map_list_including_transcript (const string& transcript);
    
    unsigned long _get_multi_map_read_count_for_combo (const string& transcript, const string& combo);

    vector<string> _get_transcripts_in_combo (string combo);

    string _create_combo (vector<string>& transcripts);

    unsigned int _get_transcript_length (const string& transcript);

    static const unsigned int MIN_EFFECTIVE_TRANS_LENGTH = 10; // large enough to keep abundance estimates from going absurdly high for short transcripts.
    

};

#endif 
