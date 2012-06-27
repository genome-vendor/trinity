#include "EM.hpp"
#include <iostream>
#include <algorithm>
#include <time.h>
#include "common.hpp"

static const string TOKEN = "^";

EM::EM (const map<string,unsigned int>& trans_seq_len_map, const unsigned int frag_length) 
    : _trans_lengths(trans_seq_len_map), _frag_length(frag_length) {

    _total_reads = 0;
    _em_round = 0;
    
}


double EM::run (unsigned int max_iterations, double min_delta_ML) {
    
    _init_theta();
    
    double prev_ML_val = 0.0;
    
    double delta = 1;

    cerr << "EM started on: " << get_num_transcripts() << " transcripts" << endl;
    
    while (true) {

        time_t start_time;
        time(&start_time);
        
        _em_round++;

        _E_step();

        _M_step();

        
        double ml = _compute_Max_Likelihood();

        if (_em_round > 1) {
            delta = ml - prev_ML_val;
        }
        
        time_t end_time;
        time(&end_time);

        time_t time_spent = end_time - start_time;

        if (DEBUG_LEVEL > 0) 
            get_and_report_results();


        
        cerr << "\rML[" << _em_round << "]: " << ml << " delta:" << delta << ", elapsed: " << time_spent << " seconds.   ";
        if (DEBUG_LEVEL > 0)
            cerr << endl;

        cerr << endl;
        
        prev_ML_val = ml;

        
        if (_em_round > 1 
            && (_em_round >= max_iterations )) { //  || delta < min_delta_ML)) {
            break;
        }
    }
    
    cerr << endl;
    
    return (prev_ML_val);
}


double EM::resume (unsigned int max_iterations, double min_delta_ML, double prev_ML_val) {
    
    // no init_theta step
    
    if (_em_round < 2) {
        stringstream errmsg;
        
        errmsg << "Error, cannot resume a run that was never started. " << stacktrace();
        throw(errmsg.str());
    }
    
    double delta;
    
    unsigned iteration_count = 0;

    while (true) {

        time_t start_time;
        time(&start_time);
        
        _em_round++;
        iteration_count++;
        
        _E_step();
        
        _M_step();
        
        double ml = _compute_Max_Likelihood();

        delta = ml - prev_ML_val;
        
        time_t end_time;
        time(&end_time);

        time_t time_spent = end_time - start_time;
        
        cerr << "\rML[" << _em_round << "]: " << ml << " delta:" << delta << ", elapsed: " << time_spent << " seconds.   ";
        
        prev_ML_val = ml;
        
        
        // get_and_report_results();
        
        if (_em_round > 1 && (iteration_count >= max_iterations  || delta < min_delta_ML)) {
            break;
        }
    }
    
        
    return (prev_ML_val);
}







void EM::_init_theta () {
    
    // ###########################
    // ## init theta:
    // ## assume multiply mapped reads equally divided among their mapped locations for init.

    vector<string> transcripts = get_all_transcripts();

    double sum_ratios = 0;

    for (vector<string>::iterator it = transcripts.begin(); it != transcripts.end(); it++) {

        string transcript = *it;
        
        double read_count_sum = 0;
        
        vector<string> combos = _get_multi_map_list_including_transcript(transcript);
        
        for (vector<string>::iterator it2 = combos.begin(); it2 != combos.end(); it2++) {
            
            string combo = *it2;
            long count = _get_multi_map_read_count_for_combo(transcript, combo);
			
            vector<string> trans_in_combo = _get_transcripts_in_combo(combo);
            int num_transcripts_share_read = trans_in_combo.size();

			//cerr << "read count for combo: " << combo << ": " << count << endl;
			
            read_count_sum += (double)count/num_transcripts_share_read;
        }

        double ratio_all_fragments = (double) read_count_sum / _total_reads;
        
        // cerr << "Initing Theta for " << transcript << ": " << ratio_all_fragments << endl;
		
        sum_ratios += ratio_all_fragments;

        _theta[transcript] = ratio_all_fragments;
    }
    
    // cerr << "SUM RATIOS: " << sum_ratios << endl;

    return;
}


vector<string> EM::get_all_transcripts() {

    vector<string> transcripts;

    for (map<string, map<string, unsigned long> >::iterator it = _trans_to_multi_map_counts.begin();
         it != _trans_to_multi_map_counts.end();
         it++) {

        string transcript = it->first;
        transcripts.push_back(transcript);
    }

    return(transcripts);
}

unsigned long EM::get_num_transcripts () {
    
    return(_trans_to_multi_map_counts.size());
}



unsigned int EM::_get_transcript_length (const string& transcript) {
    
    
    map<string,unsigned int>::const_iterator it = _trans_lengths.find(transcript);

    if (it == _trans_lengths.end()) {
        stringstream errmsg;

        errmsg << "transcript: " << transcript << " has no length value stored " << stacktrace();

        throw(errmsg.str());
    }

    return(it->second);
}


void EM::_E_step() {
  
    // cerr << "E-step()\n";
  
  
    for (map<string, map<string, unsigned long> >::iterator it = _trans_to_multi_map_counts.begin();
         it != _trans_to_multi_map_counts.end();
         it++) {

        double expected_read_count = 0;
        
        string transcript = it->first;
                
        unsigned int trans_length = _get_transcript_length(transcript);
        
        int effective_trans_length = trans_length - _frag_length + 1;
        if (effective_trans_length <= MIN_EFFECTIVE_TRANS_LENGTH) {
            effective_trans_length = MIN_EFFECTIVE_TRANS_LENGTH;
        }

        // cerr << "E: " << transcript << ", length: " << trans_length << endl;
        

        double theta_t = _get_theta(transcript);

        map<string,unsigned long> combo_map = it->second;
        
        for (map<string,unsigned long>::iterator it2 = combo_map.begin(); it2 != combo_map.end(); it2++) {
            
            string combo = it2->first;
            unsigned long combo_count = it2->second;
            
            vector<string> trans_in_combo = _get_transcripts_in_combo(combo);
            
            if (trans_in_combo.size() == 1) {
                // only current transcript uniquely mapped by reads
                expected_read_count += combo_count;
            }
            else {
                // compute partial mapping
                double numerator = ( (double) 1/effective_trans_length) * theta_t;
                //double numerator = theta_t;
                
                // denominator:  sum for all t: P(r|t) * theta(t)
                double denominator = 0.0f;
                
                for (vector<string>::iterator it3 = trans_in_combo.begin(); it3 != trans_in_combo.end(); it3++) {
                    string other_transcript = *it3;
                    
                    unsigned int other_transcript_length = _get_transcript_length(other_transcript);
                    int effective_other_transcript_length = other_transcript_length - _frag_length + 1;
                    if (effective_other_transcript_length < MIN_EFFECTIVE_TRANS_LENGTH) {
                        effective_other_transcript_length = MIN_EFFECTIVE_TRANS_LENGTH;
                    }
                    
                    double other_theta = _get_theta(other_transcript);
                    
                    double val = ( (double) 1/effective_other_transcript_length) * other_theta;
                    
                    //double val = other_theta;
                    
                    denominator += val;
    
                }
                
                

                double fractional_read_count = combo_count * (numerator/denominator);

                                
                expected_read_count += fractional_read_count;
            }
        }

        
        // cerr << "-expected read count for transcript: " << transcript << " is " << expected_read_count << endl;
        
        
        _ENt[transcript] = expected_read_count;

    }

    
    return;

}


void EM::_M_step() {

    // cerr << "M_step()" << endl;
            
    // theta = ENt(t) / sum( for all t: ENt(t) )

    for (map<string,double>::iterator it = _ENt.begin(); it != _ENt.end(); it++) { 
        
        string transcript = it->first;
        double ENt = it->second;
        
        double fragment_ratio = ENt / _total_reads;
        
        // cerr << "setting fragment ration of " << transcript << " to " << fragment_ratio << endl;
        
        _theta[transcript] = fragment_ratio;
    }

    return;
}


double EM::_compute_Max_Likelihood() {
    
    // ML = sum:t (ENt * log(theta(t)))

    double ML = 0;
    
    for (map<string,double>::iterator it = _ENt.begin(); it != _ENt.end(); it++) {
        
        string transcript = it->first;
        double ENt = it->second;
        
        double theta = _get_theta(transcript);
        
        if (theta > 0) {
            ML += ENt * log(theta);
        }
    }

    return(ML);
}

vector<t_EM_result> EM::get_results() {

    vector<string> transcripts = get_all_transcripts();
    
    vector<t_EM_result> results;

    for (vector<string>::iterator it = transcripts.begin(); it != transcripts.end(); it++) {

        string transcript = *it;

        double expected_num_frags = _get_ENt(transcript);
        double theta = _get_theta(transcript);
        unsigned int length = _get_transcript_length(transcript);
        
        t_read_counts read_counts = count_reads_mapped_to_transcript(transcript);
        
        
        double fpkm = expected_num_frags / ( (double) length/1e3) / ( (double) _total_reads/1e6);

        
        // cerr << "fpkm = " << expected_num_frags << " / (" << length << "/" << "1e3) / (" <<  _total_reads << "/1e6)" << " = " << fpkm << endl;
        
        t_EM_result result;
        result.trans_id = transcript;
        result.length = length;
        result.unique_map = read_counts.uniquely_mapped_read_count;
        result.multi_map = read_counts.multiply_mapped_read_count;
        result.expected_map = expected_num_frags;
        result.theta = theta;
        result.FPKM = fpkm;

        results.push_back(result);
    }

    return(results);
}


vector<t_EM_result>  EM::get_and_report_results () {
    
    vector<t_EM_result> results = get_results();
    
    
    cout << "#Total reads mapped: " << get_total_mapped_read_count() << endl;
    cout << "#Total transcripts examined: " << get_num_transcripts() << endl;
    cout << "#transcript\t" << "trans_length\t" << "unique_map\t" << "multi_map\t" << "Theta\t" << "EM_frag_count\t" << "FPKM" << endl;
    
    for(vector<t_EM_result>::iterator it = results.begin(); it != results.end(); it++) {
        
        t_EM_result result = *it;
        
        cout << result.trans_id << "\t"
             << result.length << "\t"
             << result.unique_map << "\t"
             << result.multi_map << "\t"
             << result.theta << "\t"
             << result.expected_map << "\t"
             << result.FPKM << endl;
    }

    return (results);
}


map<string, unsigned long>& EM::_get_combo_map_for_transcript (const string& transcript) {
    
    map<string, map<string, unsigned long> >::iterator it = _trans_to_multi_map_counts.find(transcript);
    
    if (it == _trans_to_multi_map_counts.end()) {
        stringstream errmsg;
        errmsg << "Error, transcript: " << transcript << " has no apparent read support. " << stacktrace();
        throw(errmsg.str());
    }
    
    map<string, unsigned long>& combo_map = it->second;
    
    return(combo_map);
}

map<string,  unsigned long>& EM::_get_or_create_combo_map_for_transcript (const string& transcript) {

    map<string, map<string, unsigned long> >::iterator it = _trans_to_multi_map_counts.find(transcript);
    
    if (it == _trans_to_multi_map_counts.end()) {
        
        // create one
        map<string, unsigned long> combo_map;
        _trans_to_multi_map_counts[transcript] = combo_map;

    }
        
    return(_get_combo_map_for_transcript(transcript));
}


vector<string> EM::_get_multi_map_list_including_transcript(const string& transcript) {

    map<string, unsigned long>& combo_map = _get_combo_map_for_transcript(transcript);
    
    vector<string> combos;
    
    for (map<string, unsigned long>::iterator it2 = combo_map.begin(); it2 != combo_map.end(); it2++) {

        string combo = it2->first;
        combos.push_back(combo);
    }

    return(combos);
}




unsigned long EM::_get_multi_map_read_count_for_combo(const string& transcript, const string& combo) {

 
    map<string, unsigned long>& combo_map = _get_combo_map_for_transcript(transcript);
    
    map<string, unsigned long>::iterator it = combo_map.find(combo);

    if (it == combo_map.end()) {
        stringstream errmsg;
        
        errmsg << "Error, the combo: " << combo << " is not found in multi-map for transcript: " << transcript << stacktrace();
        throw(errmsg.str());
    }

    return(it->second);
}


vector<string> EM::_get_transcripts_in_combo(string combo) {
    
    vector<string> transcripts;
    
    map<string, vector<string> >::iterator it = _combo_memoized.find(combo);
    if (it == _combo_memoized.end()) {
        // create it and store it
        
        
        string_util::tokenize(combo, transcripts, TOKEN);
        _combo_memoized[combo] = transcripts;
    }
    else {
        transcripts = it->second;
    }
    
    return(transcripts);
}


string EM::_create_combo(vector<string>& transcripts) {
    
    sort(transcripts.begin(), transcripts.end());
    
    string combo = string_util::join(transcripts, TOKEN);

    return(combo);
}


void EM::add_read(vector<string>& transcripts) {

    string combo = _create_combo(transcripts);
    
	// cerr << "Adding read: " << combo << endl;
	
    for (vector<string>::const_iterator it = transcripts.begin(); it != transcripts.end(); it++) {

        string transcript = *it;
        
        map<string, unsigned long>& combo_map = _get_or_create_combo_map_for_transcript(transcript);
        
        if (combo_map.find(combo) == combo_map.end()) {
            combo_map[combo] = 1;
        }
        else {
            combo_map[combo] = combo_map[combo] + 1;
        }
        

    }
    
    _total_reads++;

    return;
}


t_read_counts EM::count_reads_mapped_to_transcript (const string& transcript) {

    unsigned long unique_read_count = 0;
    unsigned long multi_map_count = 0;

    map<string, unsigned long>& combo_map = _get_combo_map_for_transcript(transcript);
    
    for (map<string, unsigned long>::iterator it = combo_map.begin(); it != combo_map.end(); it++) {
        
        string combo = it->first;
        unsigned long count = it->second;
        
        vector<string> transcripts = _get_transcripts_in_combo(combo);
        if (transcripts.size() == 1) {
            // uniquely mapped read
            unique_read_count += count;
        }
        else {
            // multiply mapped read
            multi_map_count += count;
        }
    }

    t_read_counts read_counts;
    read_counts.uniquely_mapped_read_count = unique_read_count;
    read_counts.multiply_mapped_read_count = multi_map_count;
    
    return(read_counts);
}


double EM::_get_theta(const string& transcript) {

    if (_theta.find(transcript) == _theta.end()) {
        
        stringstream errmsg;
        errmsg << "Error, no theta value stored for transcript: " << transcript << stacktrace();
        throw(errmsg.str());
    }
    else {
        return(_theta[transcript]);
    }
}

double EM::_get_ENt(const string& transcript) {
    
    if (_ENt.find(transcript) == _ENt.end()) {
       
        stringstream errmsg;
        errmsg << "Error, no ENt value stored for transcript: " << transcript << stacktrace();
        throw(errmsg.str());
    }
    else {
        return(_ENt[transcript]);
    }
}


unsigned long EM::get_total_mapped_read_count() {
    return(_total_reads);
}

bool EM::purge_low_percent_expressed_transcripts(float min_percent_expressed) {
    
    double sum_fpkm = 0;
    
    vector<t_EM_result> em_results = get_results();
    
    // compute sum_fpkm
    for (vector<t_EM_result>::iterator it = em_results.begin(); it != em_results.end(); it++) {
        
        t_EM_result result = *it;
        double fpkm = result.FPKM;
        sum_fpkm += fpkm;
    }

    // identify and weed out lowly expressed transcripts
    vector<string> transcripts_to_delete;
    for (vector<t_EM_result>::iterator it = em_results.begin(); it != em_results.end(); it++) {
        
        t_EM_result result = *it;
        double fpkm = result.FPKM;
        
        double percent_expressed = fpkm/sum_fpkm * 100;

        if (percent_expressed < min_percent_expressed) {
            transcripts_to_delete.push_back(result.trans_id);
        
            cerr << "-purging transcript: " << result.trans_id << " as too lowly expressed: " << percent_expressed << "%" << endl;
        }
        
    }
    
    if(! transcripts_to_delete.empty()) {
        for (vector<string>::iterator it = transcripts_to_delete.begin();
             it != transcripts_to_delete.end();
             it++) {
            
            string transcript = *it;
            
            _trans_to_multi_map_counts.erase(transcript);
            _ENt[transcript] = 0;
            _theta[transcript] = 0;
        }
        
        // reset theta values:
        _M_step();
        
        return(true);
    }
    else {
        return(false); // no transcripts to delete
    }
        
}
