#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <time.h>

#define VERSION 3.0
#define SL 35
#define MINIMIZER 32

using namespace std;

int main(int argc, char* argv[]) {
  auto start = chrono::steady_clock::now();
  srand(time(NULL));

  // Display version number.
  cout << "v." << std::fixed << std::setprecision(1) << VERSION << endl;

  string input_fasta_file;
  string output_fasta_file;

  int cf_tmp;
  opterr = 0;

  while (1) {
    static struct option long_options[] = {
        {"input-fasta-file", 1, 0, 'i'},
        {"output-fasta-file", 1, 0, 'o'},
        {0, 0, 0, 0},
    };

    int option_index = 0;
    cf_tmp = getopt_long(argc, argv, "i:o:", long_options, &option_index);

    if ((optarg != NULL) && (*optarg == '-')) {
      cf_tmp = ':';
    }

    if (cf_tmp == -1)
      break;
    switch (cf_tmp) {
    case 'i':
      input_fasta_file = optarg;
      break;
    case 'o':
      output_fasta_file = optarg;
      break;
    case ':':
      printf("Missing option for '-%s'.\n", argv[optind - 2]);
      if (long_options[option_index].has_arg == 1) {
        return 1;
      }
      break;
    case '?':
      if (optopt == 'i')
        fprintf(stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'o')
        fprintf(stderr, "Option -%c requires an argument.\n", optopt);
      else if (isprint(optopt))
        fprintf(stderr, "Unknown option `-%c'.\n", optopt);
      else
        fprintf(stderr, "Unknown option '%s'.\n", argv[optind - 1]);
      return 1;
    default:
      abort();
    }
  }

  if (argc != 5) {
    cout << "Number of supplied arguments is not correct." << endl;
    exit(0);
  }

  // Open input file.
  ifstream fin(input_fasta_file);
  if (!fin.is_open()) {
    cout << "Cannot open input FASTA file." << endl;
    exit(0);
  }
  cout << "Input file : " << input_fasta_file << endl;

  // Open output file.
  ofstream outputFile;
  outputFile.open(output_fasta_file);
  if (!outputFile.is_open()) {
    cout << "Cannot open output FASTA file." << endl;
    exit(0);
  }
  cout << "Output file : " << output_fasta_file << endl << endl;

  cout << "SL (k-mer length) = " << int(SL) << endl;
  cout << "Minimizer length = " << MINIMIZER << endl << endl;

  string line;
  string name;
  uint64_t kmer_count = 0;
  uint64_t minimizer_count = 0;

  while (std::getline(fin, line)) {
    if (line.rfind('>', 0) == 0) {
      name = line;
    }
    if (line.rfind('>', 0) != 0) {
      std::string min_str(MINIMIZER, 'Z');

      for (uint64_t i = 0; i < line.length() - MINIMIZER + 1; i++) {
        string kmer_str = line.substr(i, int(MINIMIZER));
        if (kmer_str < min_str) {
          min_str = kmer_str;
        }
      }

      if (kmer_count % 1000000 == 0) {
        cout << "-- Minimizing : " << kmer_count << " --" << endl;
      }

      // Write minimizer for a given kmer to file.
      outputFile << name << endl;
      outputFile << min_str << endl;

      minimizer_count += 1;
      kmer_count += 1;
    }
  }

  fin.close();
  outputFile.close();

  cout << endl << "Total kmers : " << kmer_count << endl;
  cout << "Total minimizers  : " << minimizer_count << endl << endl;

  // Recording end time.
  auto end = chrono::steady_clock::now();
  cout << "-- Done writing. Time so far: "
       << chrono::duration_cast<chrono::seconds>(end - start).count() << " seconds." << endl;

  return 0;
}
