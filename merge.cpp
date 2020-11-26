#include <boost/format.hpp>
#include <fstream>
#include <iostream>
#include <string>
using namespace std;
using namespace boost;


int main() {
  string lattice = "Triangular";
  string Js = "1_1_500";
  int id_min = 80000, id_max = 90000;

  string dir = (format("./out/%s/%s/") % lattice % Js).str();

  string wfile = (format("%s%s_all.txt") % dir % Js).str();
  ofstream ofs(wfile);
  for (int id = id_min; id <= id_max; id++) {
    string rfile = (format("%sslurm-%d.out") % dir % id).str();
    ifstream ifs(rfile);
    if (ifs.fail())
      continue;

    string str;
    int line = 0;
    while (getline(ifs, str)) {
      line++;
      if (line == 3)
        ofs << str << endl;
    }
  }
  return 0;
}