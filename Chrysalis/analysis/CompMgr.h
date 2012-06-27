#ifndef COMPMGR_H
#define COMPMGR_H

#include <string.h>


class ComponentFileMgr
{
 public:
  ComponentFileMgr(const string & base = ".") {
    m_base = base + "/";
    m_per = 2000;
    m_dir = m_base + "RawComps.";
  }

  string GetFileName(int index, const string & suffix) {
    int i = index / m_per;
    char tmp[1024];
    sprintf(tmp, "%s%d", m_dir.c_str(), i);
    string test = tmp;
    test += "/exists";
    FILE * p = fopen(test.c_str(), "w");
    if (p == NULL) {
      string mk = "mkdir ";
      mk += tmp;      
      // cout << "Making directory " << tmp << endl;
      system(mk.c_str());
    } else {
      fclose(p);
    }
    string out = tmp;
    sprintf(tmp, "/comp%d", index);
    out += tmp;
    out += suffix;
    return out;
  }


 private:
  string m_base;
  string m_dir;
  int m_per;

};


#endif

