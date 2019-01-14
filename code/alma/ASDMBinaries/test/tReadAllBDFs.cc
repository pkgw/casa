#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <dirent.h>

using namespace std;

// replacing boost functions not available in c++11 using stat
#include <sys/stat.h>
inline bool file_exists(const string &filename) {
  struct stat statbuf;
  return (stat(filename.c_str(),&statbuf)==0);
}

inline bool is_regular_file(const string &filename) {
  struct stat statbuf;
  bool result = (stat(filename.c_str(),&statbuf)==0);
  if (result) result = S_ISREG(statbuf.st_mode);
  return result;
}

inline bool is_directory(const string &filename) {
  struct stat statbuf;
  bool result = (stat(filename.c_str(),&statbuf)==0);
  if (result) result = S_ISDIR(statbuf.st_mode);
  return result;
}

inline off_t file_size(const string &filename) {
  struct stat statbuf;
  off_t result = 0;
  if (stat(filename.c_str(),&statbuf)==0) {
    result = statbuf.st_size;
  }
  return result;
}

int main(int argc, char* argv[])
{
  if (argc < 2)
  {
    cout << "Usage: tReadAllBDFs path\n";
    return 1;
  }

  string p (argv[1]);   // p reads clearer than argv[1] in the following code

  try
  {
    if (file_exists(p))    // does p actually exist?
    {
      if (is_regular_file(p))        // is p a regular file?
        cout << p << " size is " << file_size(p) << '\n';

      else if (is_directory(p))      // is p a directory?
      {
        cout << p << " is a directory containing:\n";

        typedef vector<string> vec;             // store paths,
        vec v;                                // so we can sort them later

	DIR *dir;
	if ((dir = opendir(p.c_str())) != NULL) {
	  struct dirent *ent;
	  while ((ent=readdir(dir)) != NULL) {
	    v.push_back(string(ent->d_name));
	  }
	}

        sort(v.begin(), v.end());             // sort, in case readdir did not proceed in sorted order

        for (vec::const_iterator it(v.begin()), it_end(v.end()); it != it_end; ++it)
        {
          cout << "   " << *it << '\n';
        }
      }
      else
        cout << p << " exists, but is neither a regular file nor a directory\n";
    }
    else
      cout << p << " does not exist\n";
  }

  catch ( exception & ex)
  {
    cout << ex.what() << endl;
  }

  return 0;
}
