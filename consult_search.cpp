#include <algorithm>           // for max
#include <bits/getopt_core.h>  // for optarg, optopt, optind, opterr
#include <bits/stdint-intn.h>  // for int64_t
#include <bits/stdint-uintn.h> // for uint64_t, uint8_t, uint32_t, uint16_t
#include <cassert>             // for assert
#include <chrono>              // for seconds, duration_cast, operator-
#include <cmath>
#include <cstdint> // for int8_t
#include <cstdio>  // for fread, fclose, fopen, fprintf, FILE
#include <cstdlib>
#include <ctype.h>
#include <ctype.h>  // for isprint
#include <dirent.h> // for closedir, dirent, opendir, readdir, DIR
#include <fstream>
#include <getopt.h> // for getopt_long, option
#include <iomanip>  // for operator<<, setprecision
#include <iostream> // for operator<<, endl, basic_ostream, ostream
#include <map>      // for map
#include <memory>   // for allocator, allocator_traits<>::value_...
#include <new>      // for bad_alloc
#include <ostream>  // for operator<<
#include <sstream>
#include <stdio.h>
#include <stdlib.h> // for exit, atoi, abort
#include <string.h> // for strcmp, strlen
#include <string>   // for string, operator+, char_traits, opera...
#include <sys/types.h>
#include <unistd.h>
#include <utility> // for pair
#include <vector>  // for vector

#define VERSION 18.0
#define KMER_LENGTH 32

#define SAVE_DISTANCES_OPT 'd'
#define MAXIMUM_DISTANCE_OPT 'D'
#define UNCLASSIFIED_OUT_OPT 0
#define CLASSIFIED_OUT_OPT 1
#define THREAD_COUNT_OPT 'T'

using namespace std;

// Prototypes.
vector<string> list_dir(const char *path);
uint8_t hd(uint64_t x, uint64_t y);
uint8_t get_enc_id(uint64_t sind, uint8_t enc_arr_id[]);
uint64_t encode_kmer_bits_reverse(const char *s, vector<int> pos);
uint64_t encode_kmer_bits(uint64_t val, vector<int8_t> shifts,
                          vector<int8_t> bits_to_grab);

void encode_kmer(const char s[], uint64_t &b_enc, uint64_t &b_sig);
void encode_kmer_reverse(const char s[], uint64_t &b_enc, uint64_t &b_sig);

void update_kmer(const char *s, uint64_t &b_enc, uint64_t &b_sig);
void update_kmer_reverse(const char *s, uint64_t &b_enc, uint64_t &b_sig);

void check_distance(uint8_t &dist, uint64_t &p, uint8_t &min_dist,
                    bool &kmer_found, bool &exact_match, uint8_t &num_matched);

void output_reads(ofstream &ofs_reads, string name, string orig_read,
                  string line_third, string line_fourth);
void output_distances(ofstream &ofs_kmer_distances, string name,
                      map<uint64_t, uint64_t> min_distances,
                      map<uint64_t, uint64_t> reverse_min_distances);

map<uint64_t, uint64_t> init_distance_map(uint64_t maximum_distance);

int main(int argc, char *argv[]) {
  auto start = chrono::steady_clock::now();

  // Display version number.
  cout << "v." << std::fixed << std::setprecision(1) << VERSION << endl;

  string input_library_dir = string();
  string query_fastq_file = string();
  string output_dir = ".";

  uint64_t k = KMER_LENGTH;
  /* Defaults. */
  uint64_t c = 1;
  uint64_t thread_count = 1;

  uint64_t maximum_distance;
  bool given_maximum_distance = false;

  bool save_distances = false;
  bool classified_out = false;
  bool unclassified_out = true;

  int cf_tmp;
  opterr = 0;

  while (1) {
    static struct option long_options[] = {
        {"input-library-dir", 1, 0, 'i'},
        {"output-dir", 1, 0, 'o'},
        {"query-fastq-file", 1, 0, 'q'},
        {"number-of-matches", 1, 0, 'c'},
        {"thread-count", 1, 0, THREAD_COUNT_OPT},
        {"save-distances", 0, 0, SAVE_DISTANCES_OPT},
        {"maximum-distance", 1, 0, MAXIMUM_DISTANCE_OPT},
        {"unclassified-out", 0, 0, UNCLASSIFIED_OUT_OPT},
        {"classified-out", 0, 0, CLASSIFIED_OUT_OPT},
        {0, 0, 0, 0},
    };

    int option_index = 0;
    cf_tmp = getopt_long(argc, argv, "i:o:q:c:", long_options, &option_index);

    if ((optarg != NULL) && (*optarg == '-')) {
      cf_tmp = ':';
    }

    if (cf_tmp == -1)
      break;
    else if (cf_tmp == SAVE_DISTANCES_OPT)
      save_distances = true;
    else if (cf_tmp == MAXIMUM_DISTANCE_OPT) {
      if (atoi(optarg) < 0) {
        cout << "Value of --maximum-distance cannot be negative." << endl;
        exit(1);
      }
      // Ignored if save_distances flag is not given.
      maximum_distance = atoi(optarg); // Default equals to p.
      given_maximum_distance = true;
    } else if (cf_tmp == CLASSIFIED_OUT_OPT)
      classified_out = true;
    else if (cf_tmp == UNCLASSIFIED_OUT_OPT)
      unclassified_out = true;
    else if (cf_tmp == THREAD_COUNT_OPT)
      thread_count = atoi(optarg); // Default is 1.
    else {
      switch (cf_tmp) {
      case 'i':
        input_library_dir = optarg;
        break;
      case 'o':
        output_dir = optarg;
        break;
      case 'q':
        query_fastq_file = optarg;
        break;
      case 'c':
        if (atoi(optarg) < 1) {
          cout << "Value of -c (--number-of-matches) cannot be smaller than 1."
               << endl;
          exit(1);
        }
        c = atoi(optarg); // Default is 1.
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
        else if (optopt == 'q')
          fprintf(stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint(optopt))
          fprintf(stderr, "Unknown option '-%c'.\n", optopt);
        else
          fprintf(stderr, "Unknown option '%s'.\n", argv[optind - 1]);
        return 1;
      default:
        abort();
      }
    }
  }

  size_t endpos = input_library_dir.find_last_not_of("/\\");
  input_library_dir = input_library_dir.substr(0, endpos + 1);

  if (argc <= 16) {
    cout << "Library directory : " << input_library_dir << endl;
    cout << "Query directory : " << query_fastq_file << endl;
    cout << "Number of threads : " << thread_count << endl;

  } else {
    cout << "Too many arguments are supplied.";
    exit(1);
  }

  if (!(classified_out || unclassified_out || save_distances)) {
    cout << "Nothing to do! Use at least one of the following flags:" << endl;
    cout << "'--classified-out' to report classified reads," << endl;
    cout << "'--unclassified-out' to report unclassified reads," << endl;
    cout << "'--save-distances' to save distances of matched k-mers." << endl;
    exit(1);
  }

  // Read map.
  string map_meta = "meta";
  string path = input_library_dir + '/' + map_meta;

  FILE *fmeta = fopen(path.c_str(), "rb");
  if (!fmeta) {
    cout << "Cannot open the metadata file of the input library!" << endl;
    return 1;
  }

  // Read parameters from input file.
  uint64_t p;
  uint64_t l;
  uint64_t h;
  uint64_t sigs_arr_size;
  uint64_t new_tag_arr_size;
  uint64_t kmer_count;
  uint64_t encli_0;
  uint64_t encli_1;
  uint64_t enc_arr_id_size;

  fread(&p, sizeof(uint64_t), 1, fmeta);
  fread(&l, sizeof(uint64_t), 1, fmeta);
  fread(&h, sizeof(uint64_t), 1, fmeta);
  fread(&sigs_arr_size, sizeof(uint64_t), 1, fmeta);
  fread(&new_tag_arr_size, sizeof(uint64_t), 1, fmeta);
  fread(&kmer_count, sizeof(uint64_t), 1, fmeta);
  fread(&encli_0, sizeof(uint64_t), 1, fmeta);
  fread(&encli_1, sizeof(uint64_t), 1, fmeta);
  fread(&enc_arr_id_size, sizeof(uint64_t), 1, fmeta);

  if (!given_maximum_distance) {
    maximum_distance = ceil(3 * p / 2) + 1;
  }
  assert(maximum_distance >= 0);

  cout << "Parameter configuration:" << endl;
  cout << "------------------------" << endl;
  cout << "k = " << k << endl;
  cout << "p = " << p << endl;
  cout << "l = " << l << endl;
  cout << "h = " << h << endl;
  cout << "c = " << c << endl;
  cout << "------------------------" << endl;
  cout << endl;

  if (save_distances) {
    cout << "Maximum distance to save is " << maximum_distance << "." << endl;
  }

  cout << endl << "k-mer statistics:" << endl;
  cout << "-----------------" << endl;
  cout << "k-mer count = " << kmer_count << endl;
  cout << "k-mer count array-0 = " << encli_0 << endl;
  cout << "k-mer count array-1 = " << encli_1 << endl;
  cout << "-----------------" << endl << endl;

  cout << "Library size information:" << endl;
  cout << "-------------------------" << endl;
  cout << "Tag array size  = " << new_tag_arr_size << endl;
  cout << "Signature array size = " << sigs_arr_size << endl;
  cout << "Encoding ID array size  = " << enc_arr_id_size << endl;
  cout << "-------------------------" << endl << endl;

  int vec_size = 0;
  vector<vector<int8_t>> shifts;
  vector<vector<int8_t>> grab_bits;

  string line;

  for (int i = 0; i < l; i++) {
    int8_t val_read;
    vector<int8_t> v;
    vector<int8_t> g;

    fread(&vec_size, sizeof(int8_t), 1, fmeta);

    for (int s = 0; s < vec_size; s++) {
      fread(&val_read, sizeof(int8_t), 1, fmeta);
      v.push_back(val_read);
    }

    shifts.push_back(v);
    v.clear();

    for (int s = 0; s < vec_size; s++) {
      fread(&val_read, sizeof(int8_t), 1, fmeta);
      g.push_back(val_read);
    }

    grab_bits.push_back(g);
    g.clear();
  }

  // Read sig rows and columns counts.
  uint64_t b;
  uint64_t sigs_row_count;
  fread(&b, sizeof(uint64_t), 1, fmeta);
  fread(&sigs_row_count, sizeof(uint64_t), 1, fmeta);

  // Read tag size, partition count, tag mask and big sig mask.
  uint64_t t;
  fread(&t, sizeof(uint64_t), 1, fmeta);
  uint64_t partitions;
  fread(&partitions, sizeof(uint64_t), 1, fmeta);
  uint64_t tag_mask;
  fread(&tag_mask, sizeof(uint64_t), 1, fmeta);
  uint64_t big_sig_mask;
  fread(&big_sig_mask, sizeof(uint64_t), 1, fmeta);

  // Display map information.
  cout << "Library information:" << endl;
  cout << "--------------------" << endl;
  cout << "Signatures row count = " << sigs_row_count << endl;
  cout << "Signatures column count = " << b * partitions << endl;
  cout << "Partitions = " << partitions << endl;
  cout << "Tag size in bits = " << t << endl;
  cout << "Tag mask = " << tag_mask << endl;
  cout << "Big sig mask = " << big_sig_mask << endl;
  cout << "--------------------" << endl << endl;

  // Read file chunk information.
  uint8_t sigf_chunks;
  fread(&sigf_chunks, sizeof(uint8_t), 1, fmeta);
  uint8_t tagf_chunks;
  fread(&tagf_chunks, sizeof(uint8_t), 1, fmeta);
  uint8_t encf_chunks;
  fread(&encf_chunks, sizeof(uint8_t), 1, fmeta);

  // Read members in file chunks.
  vector<uint64_t> sig_chunk_counts;
  vector<uint64_t> tag_chunk_counts;
  vector<uint64_t> enc_chunk_counts_0;
  vector<uint64_t> enc_chunk_counts_1;
  vector<uint64_t> enc_id_chunk_counts;

  vector<uint64_t> sig_chunk_cumcounts;
  vector<uint64_t> tag_chunk_cumcounts;
  vector<uint64_t> enc_chunk_cumcounts_0;
  vector<uint64_t> enc_chunk_cumcounts_1;
  vector<uint64_t> enc_id_chunk_cumcounts;

  uint64_t read_counts;
  uint64_t sum_read_counts = 0;

  for (int m = 0; m < sigf_chunks; m++) {
    fread(&read_counts, sizeof(uint64_t), 1, fmeta);
    sig_chunk_counts.push_back(read_counts);
    sig_chunk_cumcounts.push_back(sum_read_counts);
    sum_read_counts += read_counts;
  }

  sum_read_counts = 0;

  for (int m = 0; m < tagf_chunks; m++) {
    fread(&read_counts, sizeof(uint64_t), 1, fmeta);
    tag_chunk_counts.push_back(read_counts);
    tag_chunk_cumcounts.push_back(sum_read_counts);
    sum_read_counts += read_counts;
  }

  sum_read_counts = 0;

  for (int m = 0; m < encf_chunks; m++) {
    fread(&read_counts, sizeof(uint64_t), 1, fmeta);
    enc_chunk_counts_0.push_back(read_counts);
    enc_chunk_cumcounts_0.push_back(sum_read_counts);
    sum_read_counts += read_counts;
  }

  sum_read_counts = 0;

  for (int m = 0; m < encf_chunks; m++) {
    fread(&read_counts, sizeof(uint64_t), 1, fmeta);
    enc_chunk_counts_1.push_back(read_counts);
    enc_chunk_cumcounts_1.push_back(sum_read_counts);
    sum_read_counts += read_counts;
  }

  sum_read_counts = 0;

  for (int m = 0; m < sigf_chunks; m++) {
    fread(&read_counts, sizeof(uint64_t), 1, fmeta);
    enc_id_chunk_counts.push_back(read_counts);
    enc_id_chunk_cumcounts.push_back(sum_read_counts);
    sum_read_counts += read_counts;
  }

  fclose(fmeta);

  // Allocate signature array.
  uint32_t *sigs_arr;
  try {
    sigs_arr = new uint32_t[sigs_arr_size];
    cout << "Done memory allocation for signature array." << endl;
  } catch (bad_alloc &ba) {
    cerr << "Failed to allocate memory for signatures." << ba.what() << endl;
  }

  // Allocate tag array.
  int8_t *tag_arr;
  try {
    tag_arr = new int8_t[new_tag_arr_size];
    cout << "Done memory allocation for the tag array." << endl;
  } catch (bad_alloc &ba) {
    cerr << "Failed to allocate memory for tags." << ba.what() << endl;
  }

  // Allocate and read encoding array.
  uint64_t *encode_arr_0;
  try {
    encode_arr_0 = new uint64_t[encli_0];
    cout << "Done memory allocation for the array-0." << endl;
  } catch (bad_alloc &ba) {
    cerr << "Failed to allocate memory for array-0." << ba.what() << endl;
  }

  uint64_t *encode_arr_1;
  try {
    encode_arr_1 = new uint64_t[encli_1];
    cout << "Done memory allocation for the array-1." << endl;
  } catch (bad_alloc &ba) {
    cerr << "Failed to allocate memory for array-1." << ba.what() << endl;
  }

  // Allocate enc array id array.
  uint8_t *enc_arr_id;
  try {
    enc_arr_id = new uint8_t[enc_arr_id_size];
    cout << "Done memory allocation for the encoding array." << endl;
  } catch (bad_alloc &ba) {
    cerr << "Failed to allocate memory for encodings." << ba.what() << endl;
  }

  uint64_t total_sigs_read = 0;
  uint64_t total_tags_read = 0;
  uint64_t num_pairs_0 = 0;
  uint64_t num_pairs_1 = 0;
  uint64_t total_enc_id_read = 0;

  // Read files in parallel.
  vector<string> str_map_sig;
  vector<string> str_map_tag;
  vector<string> str_map_enc;
  vector<string> str_map_ence;
  vector<string> str_map_enc_id;

  for (int m = 0; m < sigf_chunks; m++) {
    string map_sig = "sig" + to_string(m);
    path = input_library_dir + "/" + map_sig;
    str_map_sig.push_back(path);
  }
  for (int m = 0; m < tagf_chunks; m++) {
    string map_tag = "tag" + to_string(m);
    path = input_library_dir + "/" + map_tag;
    str_map_tag.push_back(path);
  }
  for (int m = 0; m < encf_chunks; m++) {
    string map_enc = "enc" + to_string(m);
    path = input_library_dir + "/" + map_enc;
    str_map_enc.push_back(path);
  }
  for (int m = 0; m < encf_chunks; m++) {
    string map_ence = "ence" + to_string(m);
    path = input_library_dir + "/" + map_ence;
    str_map_ence.push_back(path);
  }
  for (int m = 0; m < sigf_chunks; m++) {
    string map_enc_id = "enc_id" + to_string(m);
    path = input_library_dir + "/" + map_enc_id;
    str_map_enc_id.push_back(path);
  }

  uint64_t temp_count;

// Read signature arrays.
#pragma omp parallel num_threads(thread_count) shared(sigs_arr, total_sigs_read)
  {
#pragma omp for
    for (int m = 0; m < sigf_chunks; m++) {
      FILE *f;
      f = fopen(str_map_sig[m].c_str(), "rb");
      if (!f) {
        cout << "Cannot open file for some signatures in the library directory!"
             << endl;
        exit(1);
      }

      // Read signs.
      temp_count = fread(sigs_arr + sig_chunk_cumcounts[m], sizeof(uint32_t),
                         sig_chunk_counts[m], f);
#pragma omp atomic
      total_sigs_read += temp_count;

      fclose(f);
    }
  }

// Read tag array.
#pragma omp parallel num_threads(thread_count) shared(tag_arr, total_tags_read)
  {
#pragma omp for
    for (int m = 0; m < tagf_chunks; m++) {
      FILE *ftag;
      ftag = fopen(str_map_tag[m].c_str(), "rb");
      if (!ftag) {
        cout << "Cannot open file for some tags in the library directory!"
             << endl;
        exit(1);
      }

      // Read tags.
      temp_count = fread(tag_arr + tag_chunk_cumcounts[m], sizeof(int8_t),
                         tag_chunk_counts[m], ftag);
#pragma omp atomic
      total_tags_read += temp_count;

      fclose(ftag);
    }
  }

// Read encoding array.
#pragma omp parallel num_threads(thread_count) shared(encode_arr_0, num_pairs_0)
  {
#pragma omp for
    for (int m = 0; m < encf_chunks; m++) {
      FILE *fenc;
      fenc = fopen(str_map_enc[m].c_str(), "rb");
      if (!fenc) {
        cout << "Cannot open file for some encodings in the library directory!"
             << endl;
        exit(1);
      }

      // Read encodings.
      temp_count = fread(encode_arr_0 + enc_chunk_cumcounts_0[m],
                         sizeof(uint64_t), enc_chunk_counts_0[m], fenc);
#pragma omp atomic
      num_pairs_0 += temp_count;

      fclose(fenc);
    }
  }

// Read encoding array.
#pragma omp parallel num_threads(thread_count) shared(encode_arr_1, num_pairs_1)
  {
#pragma omp for
    for (int m = 0; m < encf_chunks; m++) {
      // Open encode file.
      FILE *fence;
      fence = fopen(str_map_ence[m].c_str(), "rb");
      if (!fence) {
        cout << "Cannot open file for some encodings-e in the library!" << endl;
      }

      // Read encodings.
      temp_count = fread(encode_arr_1 + enc_chunk_cumcounts_1[m],
                         sizeof(uint64_t), enc_chunk_counts_1[m], fence);
#pragma omp atomic
      num_pairs_1 += temp_count;

      fclose(fence);
    }
  }

// Read enc_id array.
#pragma omp parallel num_threads(thread_count)                                 \
    shared(enc_arr_id, total_enc_id_read)
  {
#pragma omp for
    for (int m = 0; m < sigf_chunks; m++) {
      FILE *fenc_id;
      fenc_id = fopen(str_map_enc_id[m].c_str(), "rb");
      if (!fenc_id) {
        cout << "Cannot open file for some encodings ID in the library!"
             << endl;
        exit(1);
      }

      // Read sigs.
      temp_count = fread(enc_arr_id + enc_id_chunk_cumcounts[m],
                         sizeof(uint8_t), enc_id_chunk_counts[m], fenc_id);
#pragma omp atomic
      total_enc_id_read += temp_count;

      fclose(fenc_id);
    }
  }

  cout << endl << "Library statistics:" << endl;
  cout << "-------------------" << endl;
  cout << "Signatures read : " << total_sigs_read << endl;
  cout << "Tags read : " << total_tags_read << endl;
  cout << "Encodings ID read: " << total_enc_id_read << endl;
  cout << "Encodings read to array-0 : " << num_pairs_0 << endl;
  cout << "Encodings read to array-1 : " << num_pairs_1 << endl;
  cout << "-------------------" << endl << endl;

  auto end = chrono::steady_clock::now();
  cout << "Done reading. Now matching. Time so far: "
       << chrono::duration_cast<chrono::seconds>(end - start).count()
       << " seconds." << endl;

  /* string library_name = */
  /*     input_library_dir.substr(input_library_dir.find_last_of("/\\") + 1);
   */

  string query_fastq_truct =
      query_fastq_file.substr(query_fastq_file.find_last_of("/") + 1);
  query_fastq_truct =
      query_fastq_truct.substr(0, query_fastq_truct.find_last_of("."));
  string output_unclassified_path =
      output_dir + "/" + "unclassified-seq_" + query_fastq_truct;
  string output_classified_path =
      output_dir + "/" + "classified-seq_" + query_fastq_truct;
  string output_distances_path =
      output_dir + "/" + "kmer-distances_" + query_fastq_truct;

  // Read input fastq.
  ifstream ifs_reads_query(query_fastq_file);

  ofstream ofs_reads_unclassified;
  if (unclassified_out) {
    ofs_reads_unclassified.open(output_unclassified_path);
  }

  ofstream ofs_reads_classified;
  if (classified_out) {
    ofs_reads_classified.open(output_classified_path);
  }

  ofstream ofs_kmer_distances;
  if (save_distances) {
    ofs_kmer_distances.open(output_distances_path);
    ofs_kmer_distances << "READ_ID"
                       << "\t"
                       << "SEQ_TYPE";
    for (uint64_t i = 0; i <= maximum_distance; ++i) {
      ofs_kmer_distances << "\t" << i;
    }
    ofs_kmer_distances << endl;
  }

  uint64_t num_lines_read = 0;
  uint64_t reads_matched = 0;

#pragma omp parallel num_threads(thread_count)
  {
    while (!ifs_reads_query.eof()) {
      string name;
      string curr_read;
      string line_third;
      string line_fourth;

      uint64_t b_enc;
      uint64_t b_sig;

      uint64_t test_enc;
      uint8_t enc_arr_ind;

      int8_t tag;
      uint64_t kmer_sig;
      uint64_t big_sig_hash;
      uint64_t enc_start;
      uint64_t enc_end;

#pragma omp critical
      {
        getline(ifs_reads_query, name);
        getline(ifs_reads_query, curr_read);
        getline(ifs_reads_query, line_third);
        getline(ifs_reads_query, line_fourth);
      }

      if (ifs_reads_query.eof())
        break;

#pragma omp atomic
      num_lines_read += 4;

      uint8_t num_matched = 0;
      uint8_t num_matched_reverse = 0;

      map<uint64_t, uint64_t> min_distances =
          init_distance_map(maximum_distance);
      map<uint64_t, uint64_t> reverse_min_distances =
          init_distance_map(maximum_distance);

      string orig_read = curr_read;
      istringstream iss(curr_read);
      string token;

      while (getline(iss, token, 'N')) {
        if (token.length() >= int(k)) {
          b_enc = 0;
          b_sig = 0;

          for (uint64_t i = 0; i < token.length(); i++) {
            if (i == 0) {
              string kmer_str = token.substr(i, int(k));
              const char *ckmer = kmer_str.c_str();
              encode_kmer(ckmer, b_enc, b_sig);
              i = k - 1;

            } else {
              string kmer_str = token.substr(i, 1);
              const char *ckmer = kmer_str.c_str();
              update_kmer(ckmer, b_enc, b_sig);
            }

            bool kmer_found = false;
            bool exact_match = false;
            uint8_t min_dist = KMER_LENGTH;

            for (int64_t funci = 0; funci < l; funci++) {
              kmer_sig =
                  encode_kmer_bits(b_sig, shifts[funci], grab_bits[funci]);

              // Get first 2 bits of signature (effectively bits 28 - 27) of
              // 32 bit encoding as partition numbers.
              tag = (kmer_sig >> ((2 * h) - t)) & tag_mask;

              // Get last 26 bits of signature (effectively bits 26 -1 ) of
              // 32 bit encoding as sigs row number.
              big_sig_hash = kmer_sig & big_sig_mask;

              uint64_t f_tmp = sigs_row_count * partitions * funci;
              uint64_t s_tmp = big_sig_hash * partitions;

              if (tag_arr[f_tmp + s_tmp + tag] != -1) {
                if (tag == 0) {
                  enc_start = 0;
                } else {
                  enc_start = tag_arr[f_tmp + s_tmp + tag - 1];
                }
                enc_end = tag_arr[f_tmp + s_tmp + tag];

                for (uint64_t enc = enc_start; enc < enc_end; enc++) {
                  // Set encoding array id.
                  enc_arr_ind =
                      get_enc_id(b * f_tmp + b * s_tmp + enc, enc_arr_id);

                  if (enc_arr_ind == 0) {
                    test_enc =
                        encode_arr_0[sigs_arr[b * f_tmp + b * s_tmp + enc]];
                  } else {
                    test_enc =
                        encode_arr_1[sigs_arr[b * f_tmp + b * s_tmp + enc]];
                  }

                  uint8_t dist = hd(b_enc, test_enc);
                  check_distance(dist, p, min_dist, kmer_found, exact_match,
                                 num_matched);

                  // For each signature pointed row.
                  if (kmer_found && (!save_distances || exact_match)) {
                    break;
                  }
                }
              }
              // For each OR gate.
              if (kmer_found && (!save_distances || exact_match)) {
                break;
              }
            }
            // For each k-mer.
            if ((num_matched >= c) && (!save_distances)) {
              break;
              /* Default is equivalent to (kmer_found && save_distances). */
            } else if ((min_dist <= maximum_distance) && save_distances) {
              min_distances[min_dist]++;
            }
          }
        }
        // For each line.
        if ((num_matched >= c) && !save_distances) {
          break;
        }
      }

      // Try reverse complement.
      if ((num_matched < c) || save_distances) {
        int len = strlen(curr_read.c_str());
        char swap;

        for (int i = 0; i < len / 2; i++) {
          swap = curr_read[i];
          curr_read[i] = curr_read[len - i - 1];
          curr_read[len - i - 1] = swap;
        }

        istringstream iss(curr_read);

        while (getline(iss, token, 'N')) {
          if (token.length() >= int(k)) {
            b_enc = 0;
            b_sig = 0;

            for (uint64_t i = 0; i < token.length(); i++) {
              if (i == 0) {
                string kmer_str = token.substr(i, int(k));
                const char *ckmer = kmer_str.c_str();
                encode_kmer_reverse(ckmer, b_enc, b_sig);
                i = k - 1;
              } else {
                string kmer_str = token.substr(i, 1);
                const char *ckmer = kmer_str.c_str();
                update_kmer_reverse(ckmer, b_enc, b_sig);
              }

              bool kmer_found = false;
              bool exact_match = false;
              uint8_t min_dist = KMER_LENGTH;

              for (uint64_t funci = 0; funci < l; funci++) {
                kmer_sig =
                    encode_kmer_bits(b_sig, shifts[funci], grab_bits[funci]);

                // Get first 2 bits of signature (effectively bits 28 - 27)
                // of 32 bit encoding as partition numbers.
                tag = (kmer_sig >> ((2 * h) - t)) & tag_mask;

                // Get last 26 bits of signature (effectively bits 26 -1 )
                // of 32 bit encoding as sigs row number.
                big_sig_hash = kmer_sig & big_sig_mask;

                uint64_t f_tmp = sigs_row_count * partitions * funci;
                uint64_t s_tmp = big_sig_hash * partitions;

                if (tag_arr[f_tmp + s_tmp + tag] != -1) {
                  if (tag == 0) {
                    enc_start = 0;
                  } else {
                    enc_start = tag_arr[f_tmp + s_tmp + tag - 1];
                  }
                  enc_end = tag_arr[f_tmp + s_tmp + tag];

                  for (uint64_t enc = enc_start; enc < enc_end; enc++) {
                    enc_arr_ind =
                        get_enc_id(b * f_tmp + b * s_tmp + enc, enc_arr_id);

                    if (enc_arr_ind == 0) {
                      test_enc =
                          encode_arr_0[sigs_arr[f_tmp * b + b * s_tmp + enc]];
                    } else {
                      test_enc =
                          encode_arr_1[sigs_arr[f_tmp * b + s_tmp * b + enc]];
                    }

                    uint8_t dist = hd(b_enc, test_enc);
                    check_distance(dist, p, min_dist, kmer_found, exact_match,
                                   num_matched_reverse);

                    // For each signature pointed row.
                    if (kmer_found && (!save_distances || exact_match)) {
                      break;
                    }
                  }
                }
                // For each OR gate.
                if (kmer_found && (!save_distances || exact_match)) {
                  break;
                }
              }
              // For each k-mer.
              if ((num_matched_reverse >= c) && (!save_distances)) {
                break;
                /* Default is equivalent to (kmer_found && save_distances). */
              } else if ((min_dist <= maximum_distance) && save_distances) {
                reverse_min_distances[min_dist]++;
              }
            }
          }
          // For each line.
          if ((num_matched_reverse >= c) && !save_distances) {
            break;
          }
        }
      }

      if (save_distances) {
#pragma omp critical
        {
          output_distances(ofs_kmer_distances, name, min_distances,
                           reverse_min_distances);
        }
      }

      if ((num_matched < c) && (num_matched_reverse < c)) {
        if (unclassified_out) {
#pragma omp critical
          {
            output_reads(ofs_reads_unclassified, name, orig_read, line_third,
                         line_fourth);
          }
        }
      } else if ((num_matched >= c) || (num_matched_reverse >= c)) {
#pragma omp atomic
        reads_matched += 1;
        if (classified_out) {
#pragma omp critical
          {
            output_reads(ofs_reads_classified, name, orig_read, line_third,
                         line_fourth);
          }
        }
      }
    }
  }

  cout << query_fastq_file << " " << num_lines_read << " " << reads_matched
       << endl;

  ifs_reads_query.close();

  if (unclassified_out)
    ofs_reads_unclassified.close();
  if (classified_out)
    ofs_reads_classified.close();
  if (save_distances)
    ofs_kmer_distances.close();

  // Remember to delete array when it is done.
  delete[] sigs_arr;
  delete[] tag_arr;
  delete[] encode_arr_0;
  delete[] encode_arr_1;
  delete[] enc_arr_id;

  end = chrono::steady_clock::now();
  cout << "Done matching for all. Time so far: "
       << chrono::duration_cast<chrono::seconds>(end - start).count()
       << " seconds." << endl;

  return 0;
}

// Function definitions.

void encode_kmer(const char *s, uint64_t &b_enc, uint64_t &b_sig) {
  for (int i = 0; i < int(KMER_LENGTH); i++) {
    b_enc = b_enc << 1;
    b_sig = b_sig << 2;

    if (s[i] == 'T') {
      b_enc += 4294967297;
      b_sig += 3;
    } else if (s[i] == 'G') {
      b_enc += 4294967296;
      b_sig += 2;
    } else if (s[i] == 'C') {
      b_enc += 1;
      b_sig += 1;
    } else {
      b_enc += 0;
      b_sig += 0;
    }
  }
}

void encode_kmer_reverse(const char *s, uint64_t &b_enc, uint64_t &b_sig) {
  for (int i = 0; i < int(KMER_LENGTH); i++) {
    b_enc = b_enc << 1;
    b_sig = b_sig << 2;

    if (s[i] == 'A') {
      b_enc += 4294967297;
      b_sig += 3;
    } else if (s[i] == 'C') {
      b_enc += 4294967296;
      b_sig += 2;
    } else if (s[i] == 'G') {
      b_enc += 1;
      b_sig += 1;
    } else {
      b_enc += 0;
      b_sig += 0;
    }
  }
}

void update_kmer(const char *s, uint64_t &b_enc, uint64_t &b_sig) {
  b_enc = b_enc << 1;
  b_sig = b_sig << 2;

  // Create a mask that has set 32 bit only and set bit 32 to 0.
  uint64_t mask = 4294967297;
  b_enc = b_enc & ~mask;

  if (s[0] == 'T') {
    b_enc += 4294967297;
    b_sig += 3;
  } else if (s[0] == 'G') {
    b_enc += 4294967296;
    b_sig += 2;
  } else if (s[0] == 'C') {
    b_enc += 1;
    b_sig += 1;
  } else {
    b_enc += 0;
    b_sig += 0;
  }
}

void update_kmer_reverse(const char *s, uint64_t &b_enc, uint64_t &b_sig) {
  b_enc = b_enc << 1;
  b_sig = b_sig << 2;

  // Create a mask that has set 32 bit only and set bit 32 to 0.
  uint64_t mask = 4294967297;
  b_enc = b_enc & ~mask;

  if (s[0] == 'A') {
    b_enc += 4294967297;
    b_sig += 3;
  } else if (s[0] == 'C') {
    b_enc += 4294967296;
    b_sig += 2;
  } else if (s[0] == 'G') {
    b_enc += 1;
    b_sig += 1;
  } else {
    b_enc += 0;
    b_sig += 0;
  }
}

uint8_t hd(uint64_t x, uint64_t y) {
  uint64_t z1 = x ^ y;
  uint32_t z2 = z1 >> 32;
  uint32_t zc = z1 | z2;
  return __builtin_popcount(zc);
}

vector<string> list_dir(const char *path) {
  vector<string> userString;
  struct dirent *entry;
  DIR *dir = opendir(path);

  while ((entry = readdir(dir)) != NULL) {
    if ((strcmp(entry->d_name, "..") != 0) &&
        (strcmp(entry->d_name, ".") != 0)) {
      userString.push_back(string(path) + "/" + entry->d_name);
    }
  }

  closedir(dir);
  return (userString);
}

uint64_t encode_kmer_bits(uint64_t val, vector<int8_t> shifts,
                          vector<int8_t> bits_to_grab) {
  uint64_t res = 0;
  int i = 0;

  while (shifts[i] != -1) {
    val = val << shifts[i];
    asm("shld %b3, %2, %0"
        : "=rm"(res)
        : "0"(res), "r"(val), "ic"(bits_to_grab[i])
        : "cc");

    i++;
  }
  return uint64_t(res);
}

uint64_t encode_kmer_bits_reverse(const char *s, vector<int> pos) {
  uint64_t d = 0;
  for (int i = 0; i < pos.size(); i++) {
    d = d << 2;
    if (s[pos[i]] == 'A') {
      d += 3;
    } else if (s[pos[i]] == 'C') {
      d += 2;
    } else if (s[pos[i]] == 'G') {
      d += 1;
    } else {
      d += 0;
    }
  }
  return d;
}

uint8_t get_enc_id(uint64_t sind, uint8_t enc_arr_id[]) {
  uint64_t eind = sind >> 3;
  uint64_t ebit = sind % 8;

  uint8_t enc_arr_ind;
  enc_arr_ind = (enc_arr_id[eind] >> (7 - ebit)) & 1;

  return (enc_arr_ind);
}

void check_distance(uint8_t &dist, uint64_t &p, uint8_t &min_dist,
                    bool &kmer_found, bool &exact_match, uint8_t &num_matched) {
  if ((!kmer_found) && (dist <= p)) {
    kmer_found = true;
    num_matched += 1;
  }
  if (dist < min_dist) {
    min_dist = dist;
  }
  exact_match = dist == 0;
}

void output_reads(ofstream &ofs_reads, string name, string orig_read,
                  string line_third, string line_fourth) {
  ofs_reads << name << endl;
  ofs_reads << orig_read << endl;
  // Output separator and quality.
  ofs_reads << line_third << endl;
  ofs_reads << line_fourth << endl;
}

void output_distances(ofstream &ofs_kmer_distances, string name,
                      map<uint64_t, uint64_t> min_distances,
                      map<uint64_t, uint64_t> reverse_min_distances) {
  ofs_kmer_distances << name << "\t"
                     << "--";
  for (pair<const int, int> keyvaluepair : min_distances) {
    ofs_kmer_distances << "\t" << keyvaluepair.second;
  }
  ofs_kmer_distances << endl;

  ofs_kmer_distances << name << "\t"
                     << "rc";
  for (pair<const int, int> keyvaluepair : reverse_min_distances) {
    ofs_kmer_distances << "\t" << keyvaluepair.second;
  }
  ofs_kmer_distances << endl;
}

map<uint64_t, uint64_t> init_distance_map(uint64_t maximum_distance) {
  map<uint64_t, uint64_t> distances;
  for (uint64_t i = 0; i <= maximum_distance; ++i) {
    distances[i] = 0;
  }
  return distances;
}
