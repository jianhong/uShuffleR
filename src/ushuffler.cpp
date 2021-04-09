#include <Rcpp.h>
#include <getopt.h>
#include <string>
#include <vector>
#include <fstream>
#include <time.h>
#ifdef __cplusplus
extern "C" {
#endif
#include "ushuffle.h"
#ifdef __cplusplus
}
#endif

#ifndef USHUFFLER_CONF
#define USHUFFLER_CONF
#define USHUFFLER_VERSION "0.0.1"
#endif

using namespace Rcpp;

void print_version(){
  Rcout << "uShuffleR version: " << USHUFFLER_VERSION << std::endl;
}

void usage(){
  print_version();
  Rcout << "ushuffler [-h] [-v] [-io] [-kns]" << std::endl;
  Rcout << "where:" << std::endl;
  Rcout << "\tRequired:" << std::endl;
  Rcout << "\t\t-i,--input_file\tinput fasta file path" << std::endl;
  Rcout << "\t\t-o,--output_file\toutput fasta file path" << std::endl;
  Rcout << "\tOptional:" << std::endl;
  Rcout << "\t\t-k,--k_let\tspecifies the k-let size. Default: 2" << std::endl;
  Rcout << "\t\t-n,--number_per_seq\tspecifies the number of random sequences to generate. Default: 2" << std::endl;
  Rcout << "\t\t-s,--seed\tset the seed for random number generator" << std::endl;
  Rcout << "\t\t-h,--help\tshow this help text" << std::endl;
  Rcout << "\t\t-v,--version\tprint the version" << std::endl<< std::endl;
}

char** ou_ushuffler(char* s, int l, int k, int n){
  shuffle1(s, l, k);
  char** t = new char*[n];
  for(int i = 0; i<n; i++){
    if ((t[i] = (char*) malloc(l + 1)) == NULL) {
      stop("malloc failed\n");
    }
    t[i][l] = '\0';
    shuffle2(t[i]);
  }
  return(t);
}

// [[Rcpp::export]]
CharacterVector rushuffle(CharacterVector x, IntegerVector k, IntegerVector n) {
  int N = as<int>(n);
  int L = x.length();
  int K = as<int>(k);
  CharacterVector seq;
  for(int i=0; i<L; i++){
    std::string ins = as<std::string>(x[i]);
    char** res = ou_ushuffler((char*)ins.c_str(), ins.size(), K, N);
    for(int j=0; j<N; j++){
      seq.push_back(std::string(res[j]));
    }
  }
  return(seq);
}


int main(int argc, char **argv){
  static struct option longopts[] = {
    {"version", no_argument, 0, 'v'},
    {"help", no_argument, 0, 'h'},
    {"input_file", required_argument, 0, 'i'},
    {"output_file", required_argument, 0, 'o'},
    {"k_let", optional_argument, 0, 'k'},
    {"number_per_seq", optional_argument, 0, 'n'},
    {0, 0, 0, 0}
  };
  int opt;
  int index;
  std::string input;
  std::string output;
  int k=2;
  int n=2;

  // Retrieve the options:
  while( (opt = getopt_long(argc, argv, "i:o:k:n:vh", longopts, &index)) != -1 ){
    switch(opt){
    case 'i':
      input = optarg;
      break;
    case 'o':
      output = optarg;
      break;
    case 'k':
      k = atoi(optarg);
      break;
    case 'n':
      n = atoi(optarg);
      break;
    case 'v':
      print_version();
      break;
    case 'h':
      usage();
      break;
    case '?':
      Rcout << "Illegal option: " << opt << std::endl;
      usage();
      break;
    }
  }
  
  std::ifstream fin(input);
  if(!fin.good()){
    stop("Error in opening ", input);
  }
  std::ofstream fout(output);
  if(!fout.good()){
    stop("Error in opening ", output);
  }
  std::string line, id, seq;
  while(std::getline(fin, line)){
    if(line.empty()) continue;
    if(line[0] == '>'){
      if(!id.empty()){//create a shuffle
        char** res = ou_ushuffler((char*)seq.c_str(), seq.size(), k, n);
        for(int i=0; i<n; i++){
          fout << id << "_shuffle_" << i+1 <<std::endl;
          fout << res[i] << std::endl;
        }
      }
      id = line;
      seq.clear();
    }else{
      seq += line;
    }
  }
  
  fin.close();
  
  if(!id.empty()){//create a shuffle
    char** res = ou_ushuffler((char*)seq.c_str(), seq.size(), k, n);
    for(int i=0; i<n; i++){
      fout << id << "_shuffle_" << i+1 <<std::endl;
      fout << res[i] << std::endl;
    }
  }
  
  fout.close();
  
  return(0);
}
