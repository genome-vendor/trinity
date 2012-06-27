/* -*- mode: C++; c-basic-offset: 4; indent-tabs-mode: nil -*- */
#include <string>
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() (1)
#define omp_get_num_threads() (1)
#define omp_get_thread_num()  (0)
#endif

#include "base/CommandLineParser.h"
#include "analysis/DNAVector.h"
#include "analysis/NonRedKmerTable.h"

#include "analysis/CompMgr.h"
#include <queue>
#include <map>

#define MAX_OPEN_FILES 1000

class Assignment
{
public:
    Assignment() {
        m_comp = -1;
        m_read = -1;
    }

    void SetComp(int c) {
        m_comp = c;
    }
    void SetRead(int r) {
        m_read = r;
    }

    int Read() const {return m_read;}
    int Comp() const {return m_comp;}

    bool operator < (const Assignment & a) const {
        return (m_comp < a.m_comp);
    }

private:
    int m_comp;
    int m_read;
};



int main(int argc,char** argv)
{


    commandArg<string> aStringCmmd("-i","reads fasta");
    commandArg<string> fStringCmmd("-f","fasta input file (concatenated flat components)");
    commandArg<string> bStringCmmd("-o","output directory");
    commandArg<long> mmrCmmd("-max_mem_reads","Maximum number of reads to load into memory", -1);
    commandArg<bool> strandCmmd("-strand","strand specific data", false);
    commandArg<int> pctReadMapCmmd("-p", "percent of read kmers require mapping to component", 0);
    commandArg<bool> verboseCmmd("-verbose", "prints more status info", false);

    commandLineParser P(argc,argv);
    P.SetDescription("Assigns reads to graph components.");
    P.registerArg(aStringCmmd);
    P.registerArg(fStringCmmd);
    P.registerArg(bStringCmmd);
    P.registerArg(mmrCmmd);
    P.registerArg(strandCmmd);
    P.registerArg(pctReadMapCmmd);
    P.registerArg(verboseCmmd);
    
    P.parse();

    cerr << "-------------------------------------------" << endl
         << "---- Chrysalis: ReadsToTranscripts --------" << endl
         << "-- (Place reads on Inchworm Bundles) ------" << endl
         << "-------------------------------------------" << endl << endl;
    

    string aString = P.GetStringValueFor(aStringCmmd); // reads.fa
    string bString = P.GetStringValueFor(bStringCmmd); // out directory
    string fString = P.GetStringValueFor(fStringCmmd); // inchworm bundled.fasta
    long max_mem_reads = P.GetLongValueFor(mmrCmmd);
    bool bStrand = P.GetBoolValueFor(strandCmmd);  // strand-specific data
    int pctReadMapRequired = P.GetIntValueFor(pctReadMapCmmd); 
    
    bool VERBOSE = P.GetBoolValueFor(verboseCmmd);
    
    vecDNAVector reads, dna;

    if(max_mem_reads > 0){
        cout << "Setting maximum number of reads to load in memory to " << max_mem_reads << endl;
        reads.setMaxSeqsToRead(max_mem_reads);
    }

    cerr << "Reading bundled inchworm contigs... " << endl;
    dna.Read(fString);
    cerr << "done!" << endl;
    //test.ReverseComplement();

    int k = 25;

    ComponentFileMgr compMgr(bString);
    multimap<int,int> componentReadMap;

    unsigned long readCount = 0;
    string readCountFile = bString + "/readcounts.out";


    //-------------------------------------------------//
    //---- Build kmer index for Iworm Bundles ---------//
    //-------------------------------------------------//
    

    NonRedKmerTable kt(k);
    kt.SetUp(dna, true);
    kt.SetAllCounts(-1);

    map<int,int> component_number_mapping;

    cerr << "Assigning kmers to Iworm bundles ... ";
    #pragma omp parallel for schedule(guided,100)
    for (int i=0; i< dna.isize(); i++) {
        const DNAVector & d = dna[i];
        string bundle_acc = dna.Name(i);
        string component_no_string = bundle_acc.substr(3); // remove >s_ prefix
        int component_no = atoi(component_no_string.c_str());
        
        #pragma omp critical
        {
            component_number_mapping[i] = component_no;
            
            cerr << "\rBundle name: [" << i << "] = " << bundle_acc << ", so component_no = " << component_no << "    ";
        }


        for (int j=0; j<=d.isize()-k; j++) {
            kt.SetCount(d, j, i);  // not really counting kmers here, but instead labeling kmers according to the iworm bundles they correspond to.
        }
    }
    cerr << "done!" << endl;


    // ------------------------------------------------//
    // ------- Assign reads to Iworm Bundles ----------//
    // ------------------------------------------------//


    map<int,bool> seen_component; // open files on first write, append on later writes.

    vector<int> read_pct_mapping_info;
    read_pct_mapping_info.resize(100000, 0);
    
    cerr << "Processing reads:" << endl;
    int retVal = 0;
    unsigned long total_reads_read = 0;
    do {
        cerr << " reading another " << max_mem_reads << "... " << endl;;
        retVal = reads.Read(aString);
        cerr << "done; ";
        if (retVal > 0) {
            total_reads_read+= retVal;
            cerr << "[" << total_reads_read << "] reads analyzed for mapping." << endl;
        }
        else {
            // reached EOF
            cerr << " Finished reading reads." << endl;
        }
        
        //cerr << "read in " << abs(retVal) << " reads" << endl;
        componentReadMap.clear();
    
        #pragma omp parallel for schedule(guided,100)
        for (int i=0; i < reads.isize(); i++) {
            const DNAVector d = reads[i];
            
            // map kmers of read to Iworm bundles
            svec<int> comp;
            comp.reserve(4000);
            int num_kmer_pos = d.isize()-k + 1;
            for (int j=0; j<=d.isize()-k; j++) {
                int c = kt.GetCountReal(d, j); // the iworm bundle containing read kmer
                if (c >= 0)
                    comp.push_back(c);
            }
            if (!bStrand) {
                DNAVector dd = d;
                dd.ReverseComplement();
                for (int j=0; j<=dd.isize()-k; j++) {
                    int c = kt.GetCountReal(dd, j);
                    if (c >= 0)
                        comp.push_back(c);
                }
            }

            // find the iworm bundle with the most kmer hits
            Sort(comp);
            int best = -1;
            int max = 0;
            int run = 0;
            for (int j=1; j<comp.isize(); j++) {
                if (comp[j] != comp[j-1] || j+1 == comp.isize()) {
                    if (run > max) {
                        max = run;
                        best = comp[j-1];
                    }
                    run = 0;
                } else {
                    run++;
                }
            }
            
            //cout << "Read " << i << " maps to " << best << " with " << max << " k-mers" << endl;
            int pct_read_mapped = (int) ((float)max/num_kmer_pos*100 + 0.5);
            if(best != -1 && pct_read_mapped >= pctReadMapRequired) {
                // Note: multiple threads shouldn't try to reorder the tree
                // at the same time (which can happen during an insert)
                #pragma omp critical
                {
                    readCount++;
                    componentReadMap.insert(pair<int,int>(best,i));  // i = read_index, best = best_iworm_bundle
                    
                    // store pct read mapping info for later reporting
                    if (i > (int) read_pct_mapping_info.size()-1) {
                        read_pct_mapping_info.resize(i+100000);
                        // cerr << "resizing vec to " << i << " + 100000 " << endl;
                    }
                    read_pct_mapping_info[i] = pct_read_mapped;
                    
                }
            }
           
        } // end of read assignment to components for this round of streaming


        //-------------------------------------------------------
        // Write reads to files based on components mapped to.
        

        //cerr << "\nWriting to files... ";
        // convert sorted map to vectors for threaded processing
        vector<int> mapComponents; // list of component_ids
        vector<vector<int> > componentReads; // list of read-vectors that correspond to the component_ids above (synched by index)
        multimap<int,int>::iterator it;
        int lastComponent = -1;
        vector<int> tmpReadIDs; // temporary hold for the reads corresponding to the last component
        for(it = componentReadMap.begin(); it != componentReadMap.end(); it++){
            if(it->first != lastComponent){  // start new component, first store previous component info.
                if(tmpReadIDs.size() > 0){
                    mapComponents.push_back(lastComponent);
                    componentReads.push_back(tmpReadIDs);
                }
                lastComponent = it->first;
                tmpReadIDs.clear();
            }
            tmpReadIDs.push_back(it->second);
        }
        if(tmpReadIDs.size() > 0){ // get the last component info stored.
            mapComponents.push_back(lastComponent);
            componentReads.push_back(tmpReadIDs);
        }

        
        // write out to files in map order
        #pragma omp parallel for schedule(guided,100)
        for(int i = 0; i < (int) mapComponents.size(); i++){
            
            // each thread accesses a unique component number, so no competition between threads in writing to component-based files.
            
            int iworm_bundle_index = mapComponents[i];
            int iworm_component_no = component_number_mapping[iworm_bundle_index];

            string name = compMgr.GetFileName(iworm_component_no, ".raw.fasta");
            ofstream outFile;
            
            if(seen_component.find(iworm_component_no) == seen_component.end()) {
                // first write                
                if (VERBOSE)
                    cerr << "First write to component: " << iworm_component_no << endl;

                outFile.open(name.c_str(), ios_base::out | ios_base::trunc);
                
                #pragma omp critical
                seen_component[iworm_component_no] = true;

            } 
            else {
                // subsequent write as append
                if (VERBOSE)
                    cerr << "Second write as append to component: " << iworm_component_no << endl;
                
                outFile.open(name.c_str(), ios_base::out | ios_base::app);                            
            }
            for(unsigned int j = 0; j < componentReads[i].size(); j++){
                int read_index = componentReads[i][j];
                outFile << reads.Name(read_index) << " " << read_pct_mapping_info[read_index] << "%" << endl;
                outFile << reads[read_index].AsString() << endl;
            }
            outFile.close();
        }
        
        //cerr << "done\n";
        //clear out read component mappings
    } while (retVal > 0);
    
    cerr << "\nDone" << endl; // don't overwrite last read count

    FILE * pReadCount = fopen(readCountFile.c_str(), "w");
    fprintf(pReadCount, "%lu\n", readCount/2);
    fclose(pReadCount);

    return 0;
}

