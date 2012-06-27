#include "analysis/NonRedKmerTable.h"

bool Regular(char c) {
  if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
    return true;
  else
    return false;
}


// Restricts all k-mers to what's in here.
void NonRedKmerTable::SetUp(const vecDNAVector & templ, bool noNs)
{
  int i, j;
  
  // cerr << "vecDNAVector length: " << templ.isize() << endl;
  
  int l = 0;
  for (i=0; i<templ.isize(); i++)
    l += templ[i].isize()-m_k+1;
  
  // cerr << "Allocating: " << l << endl;
  m_data.resize(l);

  l = 0;
  for (i=0; i<templ.isize(); i++) {
    const DNAVector & d = templ[i];
    for (j=0; j<=d.isize()-m_k; j++) {      
      d.Substring(m_data[l], j, m_k);
      
      if (noNs) {
	bool bGood = true;
	for (int x=0; x<m_k; x++) {
	  if (!Regular((m_data[l])[x]))
	    bGood = false;
	}
	if (!bGood)
	  continue;
      }
      


      //cout << m_data[l] << endl;
      l++;
    }
  }
  // cerr << "Resizing: " << l << endl;
  m_data.resize(l);

  UniqueSort(m_data);
  m_counts.resize(m_data.isize(), 0);
  // cerr << "Done." << endl;
  
 
}

void NonRedKmerTable::AddData(const vecDNAVector & data)
{
  int i, j;

  cerr << "NonRedKmerTable::AddData() - adding Reads: " << data.isize() << endl;
  for (i=0; i<data.isize(); i++) {
    
    if (i % 1000 == 0) {
      cerr << 100.*(double)i/(double)data.isize() << " %" << endl;
    }
    const DNAVector & d = data[i];
    for (j=0; j<=d.isize()-m_k; j++) {      
      string s; 
      d.Substring(s, j, m_k);
	

      int k = BinSearch(m_data, s);
      if (k < 0)
	continue;
      m_counts[k]++;
    }
  }

}

void NonRedKmerTable::AddData(vecDNAVectorStream & data)
{
    int i, j;
    
    long count_reads_read = 0;
    
    cerr << "Reads: (unknown count: streaming-mode)" << endl;
    while (true) {
        const DNAVector & d = data.Next();
        
        if (d.isize() == 0)
            break;
        
        //cerr << "read another: " << d.isize() << " reads." << endl;
        
        //count_reads_read += d.isize();
        count_reads_read++;
        
        if (count_reads_read % 1000 == 0)
            cerr << "\r[" << count_reads_read << "] reads read.  ";
        
        for (j=0; j<=d.isize()-m_k; j++) {      
            string s; 
            d.Substring(s, j, m_k);
            
            
            int k = BinSearch(m_data, s);
            // cerr << "Looking up kmer: " << s << " found " << k << " times." << endl;
            if (k < 0)
                continue;
            
            m_counts[k]++;
        }
    }
    cerr << endl;
    
    
}

