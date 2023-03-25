#include <algorithm>
#include <bits/getopt_core.h>
#include <bits/stdint-intn.h>
#include <bits/stdint-uintn.h>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <ctype.h>
#include <dirent.h>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <new>
#include <ostream>
#include <random>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <unordered_map>
#include <utility>
#include <vector>

#define KMER_LENGTH 32
#define TAX_ID_LIMIT 65534

#define SAVE_DISTANCES_OPT 'd'
#define MAXIMUM_DISTANCE_OPT 'D'
#define UNCLASSIFIED_OUT_OPT 0
#define CLASSIFIED_OUT_OPT 1
#define THREAD_COUNT_OPT 'T'
#define SAVE_MATCHES_OPT 'M'
#define INIT_ID_OPT 'I'
#define UPDATE_ID_OPT 'U'
#define TAXONOMY_LOOKUP_PATH_OPT 'A'
#define FILENAME_MAP_PATH_OPT 'F'

#define READ_PARALLELISM
/* #define FILE_PARALLELISM */

using namespace std;

// Prototypes.
vector<string> list_dir(const char *path);
uint8_t hd(uint64_t x, uint64_t y);
uint8_t get_encid(uint64_t sind, uint8_t encid_arr[]);
uint64_t encode_kmer_bits_reverse(const char *s, vector<int> pos);
uint64_t encode_kmer_bits(uint64_t val, vector<int8_t> shifts, vector<int8_t> bits_to_grab);

void encode_kmer(const char s[], uint64_t &b_enc, uint64_t &b_sig);
void encode_kmer_reverse(const char s[], uint64_t &b_enc, uint64_t &b_sig);

void update_kmer(const char *s, uint64_t &b_enc, uint64_t &b_sig);
void update_kmer_reverse(const char *s, uint64_t &b_enc, uint64_t &b_sig);

void check_distance(uint8_t &dist, uint64_t &p, uint8_t &min_dist, bool &kmer_found, bool &closest_match,
                    bool &exact_match, uint8_t &num_matched);

void output_reads(ofstream &ofs_reads, string name, string orig_read, string line_third, string line_fourth);
void output_distances(ofstream &ofs_kmer_distances, string name, map<uint64_t, uint64_t> min_distances,
                      map<uint64_t, uint64_t> reverse_min_distances);
void output_matches(ofstream &ofs_match_information, string name, vector<uint16_t> match_cIDs,
                    vector<uint16_t> match_distances, vector<uint16_t> rc_match_taxIDs,
                    vector<uint16_t> rc_match_distances, map<uint16_t, uint64_t> &cID_to_taxID);

void read_filename_map(string filename, map<string, uint64_t> &filename_to_taxID);
void read_taxonomy_lookup(string filepath, unordered_map<uint16_t, vector<uint16_t>> &taxonomy_lookup,
                          map<uint64_t, uint16_t> &taxID_to_cID, map<uint16_t, uint64_t> &cID_to_taxID);

void update_kmer_cID(uint16_t cID_arr_0[], uint16_t cID_arr_1[], uint16_t count_arr_0[],
                     uint16_t count_arr_1[], vector<bool> &seen_0, vector<bool> &seen_1, uint8_t encid,
                     uint32_t encoding_ix, uint16_t filename_cID,
                     unordered_map<uint16_t, vector<uint16_t>> &taxonomy_lookup);
void update_kmer_count(uint16_t count_arr_0[], uint16_t count_arr_1[], vector<bool> &seen_0,
                       vector<bool> &seen_1, uint8_t encid, uint32_t encoding_ix);

map<uint64_t, uint64_t> init_distance_map(uint64_t maximum_distance);

int main(int argc, char *argv[]) {
  /* omp_set_nested(1); // Enable nested parallelism */
  /* omp_set_max_active_levels(2); */
  auto start = chrono::steady_clock::now();

  string input_library_dir = string();
  string output_result_dir = ".";
  string taxonomy_lookup_path;
  string filename_map_path;
  char *query_path = NULL;

  uint64_t k = KMER_LENGTH;
  /* Defaults. */
  uint64_t c = 1;
  uint64_t thread_count = 1;

  uint64_t maximum_distance;
  bool given_maximum_distance = false;

  bool save_distances = false;
  bool classified_out = false;
  bool unclassified_out = false;
  bool update_ID = false;
  bool init_ID = false;
  bool save_matches = false;

  int cf_tmp;
  opterr = 0;

  while (1) {
    static struct option long_options[] = {
        {"input-library-dir", 1, 0, 'i'},
        {"output-result-dir", 1, 0, 'o'},
        {"query-path", 1, 0, 'q'},
        {"number-of-matches", 1, 0, 'c'},
        {"thread-count", 1, 0, THREAD_COUNT_OPT},
        {"save-distances", 0, 0, SAVE_DISTANCES_OPT},
        {"maximum-distance", 1, 0, MAXIMUM_DISTANCE_OPT},
        {"unclassified-out", 0, 0, UNCLASSIFIED_OUT_OPT},
        {"classified-out", 0, 0, CLASSIFIED_OUT_OPT},
        {"init-ID", 0, 0, INIT_ID_OPT},
        {"update-ID", 0, 0, UPDATE_ID_OPT},
        {"save-matches", 0, 0, SAVE_MATCHES_OPT},
        {"taxonomy-lookup-path", 1, 0, TAXONOMY_LOOKUP_PATH_OPT},
        {"filename-map-path", 1, 0, FILENAME_MAP_PATH_OPT},
        {0, 0, 0, 0},
    };

    int option_ID = 0;
    cf_tmp = getopt_long(argc, argv, "i:o:q:c:", long_options, &option_ID);

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
    else if (cf_tmp == UPDATE_ID_OPT)
      update_ID = true;
    else if (cf_tmp == INIT_ID_OPT)
      init_ID = true;
    else if (cf_tmp == SAVE_MATCHES_OPT)
      save_matches = true;
    else if (cf_tmp == THREAD_COUNT_OPT)
      thread_count = atoi(optarg); // Default is 1.
    else if (cf_tmp == TAXONOMY_LOOKUP_PATH_OPT)
      taxonomy_lookup_path = optarg;
    else if (cf_tmp == FILENAME_MAP_PATH_OPT)
      filename_map_path = optarg;
    else {
      switch (cf_tmp) {
      case 'i':
        input_library_dir = optarg;
        break;
      case 'o':
        output_result_dir = optarg;
        break;
      case 'q':
        // Can be a directory path or a file path.
        query_path = optarg;
        break;
      case 'c':
        if (atoi(optarg) < 1) {
          cout << "Value of -c (--number-of-matches) cannot be smaller than 1." << endl;
          exit(1);
        }
        c = atoi(optarg); // Default is 1.
        break;
      case ':':
        printf("Missing option for '-%s'.\n", argv[optind - 2]);
        if (long_options[option_ID].has_arg == 1) {
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
  if ((init_ID && update_ID) || (init_ID && save_matches) || (update_ID && save_matches)) {
    cout << "Can only do one of initialize, update or save at a time." << endl;
    exit(1);
  }

  size_t endpos = input_library_dir.find_last_not_of("/\\");
  input_library_dir = input_library_dir.substr(0, endpos + 1);

  if (argc <= 20) {
    cout << "Input library directory : " << input_library_dir << endl;
    cout << "Output result directory : " << output_result_dir << endl;
    cout << "Query path : " << query_path << endl;
    cout << "Number of threads : " << thread_count << endl;

  } else {
    cout << "Too many arguments are supplied.";
    exit(1);
  }

  if (!(classified_out || unclassified_out || save_distances || init_ID || update_ID || save_matches)) {
    cout << "Nothing to do! Use at least one of the following flags:" << endl;
    cout << "'--classified-out' to report classified reads," << endl;
    cout << "'--unclassified-out' to report unclassified reads," << endl;
    cout << "'--init-ID' to initialize taxonomic ID array." << endl;
    cout << "'--update-ID' to compute taxonomic labels for k-mers." << endl;
    cout << "'--save-matches' to save matched labels and distances." << endl;
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
  float alpha;
  uint64_t h;
  uint64_t sigs_arr_size;
  uint64_t new_tag_arr_size;
  uint64_t kmer_count;
  uint64_t encli_0;
  uint64_t encli_1;
  uint64_t encid_arr_size;

  fread(&p, sizeof(uint64_t), 1, fmeta);
  fread(&l, sizeof(uint64_t), 1, fmeta);
  fread(&alpha, sizeof(float), 1, fmeta);
  fread(&h, sizeof(uint64_t), 1, fmeta);

  fread(&sigs_arr_size, sizeof(uint64_t), 1, fmeta);
  fread(&new_tag_arr_size, sizeof(uint64_t), 1, fmeta);
  fread(&kmer_count, sizeof(uint64_t), 1, fmeta);
  fread(&encli_0, sizeof(uint64_t), 1, fmeta);
  fread(&encli_1, sizeof(uint64_t), 1, fmeta);
  fread(&encid_arr_size, sizeof(uint64_t), 1, fmeta);

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

  cout << endl
       << "k-mer statistics:" << endl;
  cout << "-----------------" << endl;
  cout << "k-mer count = " << kmer_count << endl;
  cout << "k-mer count array-0 = " << encli_0 << endl;
  cout << "k-mer count array-1 = " << encli_1 << endl;
  cout << "-----------------" << endl
       << endl;

  cout << "Library size information:" << endl;
  cout << "-------------------------" << endl;
  cout << "Tag array size  = " << new_tag_arr_size << endl;
  cout << "Signature array size = " << sigs_arr_size << endl;
  cout << "Encoding ID array size  = " << encid_arr_size << endl;
  cout << "-------------------------" << endl
       << endl;

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
  cout << "--------------------" << endl
       << endl;

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
  vector<uint64_t> encid_chunk_counts;

  vector<uint64_t> sig_chunk_cumcounts;
  vector<uint64_t> tag_chunk_cumcounts;
  vector<uint64_t> enc_chunk_cumcounts_0;
  vector<uint64_t> enc_chunk_cumcounts_1;
  vector<uint64_t> encid_chunk_cumcounts;

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
    encid_chunk_counts.push_back(read_counts);
    encid_chunk_cumcounts.push_back(sum_read_counts);
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
  uint64_t *encode_arr_1;
  try {
    encode_arr_0 = new uint64_t[encli_0];
    cout << "Done memory allocation for the array-0." << endl;
  } catch (bad_alloc &ba) {
    cerr << "Failed to allocate memory for array-0." << ba.what() << endl;
  }
  try {
    encode_arr_1 = new uint64_t[encli_1];
    cout << "Done memory allocation for the array-1." << endl;
  } catch (bad_alloc &ba) {
    cerr << "Failed to allocate memory for array-1." << ba.what() << endl;
  }

  uint16_t *cID_arr_0;
  uint16_t *cID_arr_1;
  if (init_ID || update_ID || save_matches) {
    try {
      cID_arr_0 = new uint16_t[encli_0];
      cout << "Done memory allocation for the ID array-0." << endl;
    } catch (bad_alloc &ba) {
      cerr << "Failed to allocate memory for ID array-0." << ba.what() << endl;
    }
    try {
      cID_arr_1 = new uint16_t[encli_1];
      cout << "Done memory allocation for the ID array-1." << endl;
    } catch (bad_alloc &ba) {
      cerr << "Failed to allocate memory for ID array-1." << ba.what() << endl;
    }
  }

  uint16_t *count_arr_0;
  uint16_t *count_arr_1;
  if (init_ID || update_ID) {
    try {
      count_arr_0 = new uint16_t[encli_0];
      cout << "Done memory allocation for the count array-1." << endl;
    } catch (bad_alloc &ba) {
      cerr << "Failed to allocate memory for count array-1." << ba.what() << endl;
    }
    try {
      count_arr_1 = new uint16_t[encli_1];
      cout << "Done memory allocation for the count array-1." << endl;
    } catch (bad_alloc &ba) {
      cerr << "Failed to allocate memory for count array-1." << ba.what() << endl;
    }
  }

  // Allocate enc array id array.
  uint8_t *encid_arr;
  try {
    encid_arr = new uint8_t[encid_arr_size];
    cout << "Done memory allocation for the encoding array." << endl;
  } catch (bad_alloc &ba) {
    cerr << "Failed to allocate memory for encodings." << ba.what() << endl;
  }

  cout << "Memory allocation is completed." << endl
       << endl;

  struct stat s_query_path;
  vector<string> query_file_list;
  bool is_query_dir = false;
  if (stat(query_path, &s_query_path) == 0) {
    if (s_query_path.st_mode & S_IFDIR) {
      query_file_list = list_dir(query_path);
      is_query_dir = true;
    } else if (s_query_path.st_mode & S_IFREG) {
      query_file_list.push_back(query_path);
    } else {
      cout << "Filetype in the given query path is not recognized." << endl;
      exit(1);
    }
  } else {
    cout << "Given query path is not valid!" << endl;
    exit(1);
  }

  uint64_t total_sigs_read = 0;
  uint64_t total_tags_read = 0;
  uint64_t num_pairs_0 = 0;
  uint64_t num_pairs_1 = 0;
  uint64_t total_encid_read = 0;

  // Read files in parallel.
  vector<string> str_map_sig;
  vector<string> str_map_tag;
  vector<string> str_map_cID;
  vector<string> str_map_enc;
  vector<string> str_map_count;
  vector<string> str_map_cIDe;
  vector<string> str_map_ence;
  vector<string> str_map_counte;
  vector<string> str_map_encid;

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
    string map_cID = "cID" + to_string(m);
    path = input_library_dir + "/" + map_cID;
    str_map_cID.push_back(path);
    string map_enc = "enc" + to_string(m);
    path = input_library_dir + "/" + map_enc;
    str_map_enc.push_back(path);
    string map_count = "count" + to_string(m);
    path = input_library_dir + "/" + map_count;
    str_map_count.push_back(path);

    string map_cIDe = "cIDe" + to_string(m);
    path = input_library_dir + "/" + map_cIDe;
    str_map_cIDe.push_back(path);
    string map_ence = "ence" + to_string(m);
    path = input_library_dir + "/" + map_ence;
    str_map_ence.push_back(path);
    string map_counte = "counte" + to_string(m);
    path = input_library_dir + "/" + map_counte;
    str_map_counte.push_back(path);
  }
  for (int m = 0; m < sigf_chunks; m++) {
    string map_encid = "encid" + to_string(m);
    path = input_library_dir + "/" + map_encid;
    str_map_encid.push_back(path);
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
        cout << "Cannot open file for signatures in the library directory!" << endl;
        exit(1);
      }

      // Read signs.
      temp_count = fread(sigs_arr + sig_chunk_cumcounts[m], sizeof(uint32_t), sig_chunk_counts[m], f);
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
        cout << "Cannot open file for tags in the library directory!" << endl;
        exit(1);
      }

      // Read tags.
      temp_count = fread(tag_arr + tag_chunk_cumcounts[m], sizeof(int8_t), tag_chunk_counts[m], ftag);
#pragma omp atomic
      total_tags_read += temp_count;

      fclose(ftag);
    }
  }

// Read encoding array-0.
#pragma omp parallel num_threads(thread_count) shared(encode_arr_0, num_pairs_0)
  {
#pragma omp for
    for (int m = 0; m < encf_chunks; m++) {
      FILE *fenc;
      fenc = fopen(str_map_enc[m].c_str(), "rb");
      if (!fenc) {
        cout << "Cannot open file for encodings in the library directory!" << endl;
        exit(1);
      }

      // Read encodings.
      temp_count =
          fread(encode_arr_0 + enc_chunk_cumcounts_0[m], sizeof(uint64_t), enc_chunk_counts_0[m], fenc);
#pragma omp atomic
      num_pairs_0 += temp_count;

      fclose(fenc);
    }
  }
// Read encoding array-1.
#pragma omp parallel num_threads(thread_count) shared(encode_arr_1, num_pairs_1)
  {
#pragma omp for
    for (int m = 0; m < encf_chunks; m++) {
      // Open encode file.
      FILE *fence;
      fence = fopen(str_map_ence[m].c_str(), "rb");
      if (!fence) {
        cout << "Cannot open file for encodings in the library!" << endl;
      }

      // Read encodings.
      temp_count =
          fread(encode_arr_1 + enc_chunk_cumcounts_1[m], sizeof(uint64_t), enc_chunk_counts_1[m], fence);
#pragma omp atomic
      num_pairs_1 += temp_count;

      fclose(fence);
    }
  }

  if (update_ID || save_matches) {
    // Read cIDs-0.
#pragma omp parallel num_threads(thread_count) shared(cID_arr_0)
    {
#pragma omp for
      for (int m = 0; m < encf_chunks; m++) {
        FILE *fcID;
        fcID = fopen(str_map_cID[m].c_str(), "rb");
        if (!fcID) {
          cout << "Cannot open file for k-mer cIDs in the library directory!" << endl;
          exit(1);
        }

        // Read k-mer cIDs.
        fread(cID_arr_0 + enc_chunk_cumcounts_0[m], sizeof(uint16_t), enc_chunk_counts_0[m], fcID);
        fclose(fcID);
      }
    }
    // Read cIDs-1.
#pragma omp parallel num_threads(thread_count) shared(cID_arr_1)
    {
#pragma omp for
      for (int m = 0; m < encf_chunks; m++) {
        FILE *fcID;
        fcID = fopen(str_map_cIDe[m].c_str(), "rb");
        if (!fcID) {
          cout << "Cannot open file for k-mer cIDs in the library directory!" << endl;
          exit(1);
        }

        // Read k-mer cIDs.
        fread(cID_arr_1 + enc_chunk_cumcounts_1[m], sizeof(uint16_t), enc_chunk_counts_1[m], fcID);
        fclose(fcID);
      }
    }
  }

  if (update_ID) {
    // Read counts-0.
#pragma omp parallel num_threads(thread_count) shared(count_arr_0)
    {
#pragma omp for
      for (int m = 0; m < encf_chunks; m++) {
        FILE *fcount;
        fcount = fopen(str_map_count[m].c_str(), "rb");
        if (!fcount) {
          cout << "Cannot open file for count in the library directory!" << endl;
          exit(1);
        }

        // Read k-mer counts.
        fread(count_arr_0 + enc_chunk_cumcounts_0[m], sizeof(uint16_t), enc_chunk_counts_0[m], fcount);
        fclose(fcount);
      }
    }
    // Read counts-1.
#pragma omp parallel num_threads(thread_count) shared(count_arr_1)
    {
#pragma omp for
      for (int m = 0; m < encf_chunks; m++) {
        FILE *fcounte;
        fcounte = fopen(str_map_counte[m].c_str(), "rb");
        if (!fcounte) {
          cout << "Cannot open file for counts in the library directory!" << endl;
          exit(1);
        }

        // Read k-mer counts.
        fread(count_arr_1 + enc_chunk_cumcounts_1[m], sizeof(uint16_t), enc_chunk_counts_1[m], fcounte);
        fclose(fcounte);
      }
    }
  }

  // Initialize class IDs and counts to 0.
  if (init_ID) {
    for (uint64_t j = 0; j < encli_0; j++) {
      cID_arr_0[j] = 0;
      count_arr_0[j] = 0;
    }
    for (uint64_t j = 0; j < encli_1; j++) {
      cID_arr_1[j] = 0;
      count_arr_1[j] = 0;
    }
  }

  map<string, uint64_t> filename_to_taxID;
  map<string, uint16_t> filename_to_cID;
  unordered_map<uint16_t, vector<uint16_t>> taxonomy_lookup;
  map<uint64_t, uint16_t> taxID_to_cID;
  map<uint16_t, uint64_t> cID_to_taxID;
  if (update_ID || save_matches) {
    // Read filename look-up table taxonomy IDs.
    read_filename_map(filename_map_path, filename_to_taxID);
    // Read taxonomy lookup table for computing common ancestors.
    read_taxonomy_lookup(taxonomy_lookup_path, taxonomy_lookup, taxID_to_cID, cID_to_taxID);
    for (auto const &kv : filename_to_taxID)
      filename_to_cID[kv.first] = taxID_to_cID[kv.second];
  }

// Read encid array.
#pragma omp parallel num_threads(thread_count) shared(encid_arr, total_encid_read)
  {
#pragma omp for
    for (int m = 0; m < sigf_chunks; m++) {
      FILE *fencid;
      fencid = fopen(str_map_encid[m].c_str(), "rb");
      if (!fencid) {
        cout << "Cannot open file for encoding IDs in the library!" << endl;
        exit(1);
      }

      // Read sigs.
      temp_count = fread(encid_arr + encid_chunk_cumcounts[m], sizeof(uint8_t), encid_chunk_counts[m], fencid);
#pragma omp atomic
      total_encid_read += temp_count;

      fclose(fencid);
    }
  }

  cout << "Library statistics:" << endl;
  cout << "-------------------" << endl;
  cout << "Signatures read : " << total_sigs_read << endl;
  cout << "Tags read : " << total_tags_read << endl;
  cout << "Encoding IDs read: " << total_encid_read << endl;
  cout << "Encodings read to array-0 : " << num_pairs_0 << endl;
  cout << "Encodings read to array-1 : " << num_pairs_1 << endl;
  cout << "-------------------" << endl
       << endl;

  auto end = chrono::steady_clock::now();
  cout << "Done reading. Now matching. Time so far: "
       << chrono::duration_cast<chrono::seconds>(end - start).count() << " seconds." << endl;

  int counter_files = 1;

#ifdef FILE_PARALLELISM
#pragma omp parallel for schedule(dynamic) num_threads(thread_count)
#endif
  for (int fc_ix = 0; fc_ix < query_file_list.size(); ++fc_ix) {
    string query_file_path;
#ifdef FILE_PARALLELISM
#pragma omp critical
#endif
    {
      query_file_path = query_file_list[fc_ix];
      cout << counter_files << "/" << query_file_list.size() << endl;
      counter_files++;
    }

    vector<bool> seen_0;
    vector<bool> seen_1;

    if (init_ID || update_ID) {
      seen_0.resize(encli_0, false);
      seen_1.resize(encli_1, false);
    }

    string query_fastq_truct = query_file_path.substr(query_file_path.find_last_of("/") + 1);
    query_fastq_truct = query_fastq_truct.substr(0, query_fastq_truct.find_last_of("."));
    string output_unclassified_path = output_result_dir + "/" + "unclassified-seq_" + query_fastq_truct;
    string output_classified_path = output_result_dir + "/" + "classified-seq_" + query_fastq_truct;
    string output_distances_path = output_result_dir + "/" + "kmer-distances_" + query_fastq_truct;
    string output_matches_path = output_result_dir + "/" + "match-info_" + query_fastq_truct;

    uint16_t filename_cID = filename_to_cID[query_fastq_truct];

    // Read input fastq.
    ifstream ifs_reads_query(query_file_path);

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

    ofstream ofs_match_information;
    if (save_matches) {
      ofs_match_information.open(output_matches_path);
    }

    uint64_t num_lines_read = 0;
    uint64_t reads_matched = 0;

#ifdef READ_PARALLELISM
#pragma omp parallel num_threads(thread_count)
#endif
    {
      while (ifs_reads_query) {
        string name;
        string curr_read;
        string line_third;
        string line_fourth;

        uint64_t b_enc;
        uint64_t b_sig;

        uint64_t test_enc;
        uint8_t encid;

        int8_t tag;
        uint64_t kmer_sig;
        uint64_t big_sig_hash;
        uint64_t enc_start;
        uint64_t enc_end;

#ifdef READ_PARALLELISM
#pragma omp critical
#endif
        {
          getline(ifs_reads_query, name);
          getline(ifs_reads_query, curr_read);
          getline(ifs_reads_query, line_third);
          getline(ifs_reads_query, line_fourth);
        }

        if (ifs_reads_query.eof())
          break;

#ifdef READ_PARALLELISM
#pragma omp atomic
#endif
        num_lines_read += 4;

        uint8_t num_matched = 0;
        uint8_t num_matched_reverse = 0;

        map<uint64_t, uint64_t> min_distances = init_distance_map(maximum_distance);
        map<uint64_t, uint64_t> reverse_min_distances = init_distance_map(maximum_distance);

        vector<uint16_t> match_distances;
        vector<uint16_t> match_cIDs;
        vector<uint16_t> rc_match_distances;
        vector<uint16_t> rc_match_taxIDs;

        string orig_read = curr_read;
        istringstream iss(curr_read);
        string token;

        while (getline(iss, token, 'N')) {
          if (token.length() >= int(k)) {
            b_enc = 0;
            b_sig = 0;

            for (uint64_t i = 0; i < token.length(); i++) {
              string kmer_str;
              if (i == 0) {
                kmer_str = token.substr(i, int(k));
                const char *ckmer = kmer_str.c_str();
                encode_kmer(ckmer, b_enc, b_sig);
                i = k - 1;

              } else {
                kmer_str = token.substr(i, 1);
                const char *ckmer = kmer_str.c_str();
                update_kmer(ckmer, b_enc, b_sig);
              }

              bool kmer_found = false;
              bool closest_match = false;
              bool exact_match = false;
              uint8_t min_dist = KMER_LENGTH;
              uint16_t kmer_ID;

              for (int64_t funci = 0; funci < l; funci++) {
                kmer_sig = encode_kmer_bits(b_sig, shifts[funci], grab_bits[funci]);

                tag = (kmer_sig >> ((2 * h) - t)) & tag_mask;
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
                    encid = get_encid(b * f_tmp + b * s_tmp + enc, encid_arr);

                    uint32_t encoding_ix = sigs_arr[b * f_tmp + b * s_tmp + enc];
                    if (encid == 0) {
                      test_enc = encode_arr_0[encoding_ix];
                    } else {
                      test_enc = encode_arr_1[encoding_ix];
                    }

                    uint8_t dist = hd(b_enc, test_enc);
                    check_distance(dist, p, min_dist, kmer_found, closest_match, exact_match, num_matched);

                    if (update_ID && exact_match) {
#pragma omp critical
                      {
                        update_kmer_cID(cID_arr_0, cID_arr_1, count_arr_0, count_arr_1, seen_0, seen_1, encid,
                                        encoding_ix, filename_cID, taxonomy_lookup);
                      }
                    }
                    if (init_ID && exact_match) {
#pragma omp critical
                      {
                        update_kmer_count(count_arr_0, count_arr_1, seen_0, seen_1, encid, encoding_ix);
                      }
                    }
                    if (closest_match && save_matches) {
                      if (encid == 0) {
                        kmer_ID = cID_arr_0[encoding_ix];
                      } else {
                        kmer_ID = cID_arr_1[encoding_ix];
                      }
                    }
                    // For each signature pointed row.
                    if ((!init_ID) && (!update_ID) && kmer_found &&
                        (exact_match || (!save_distances && !save_matches))) {
                      break;
                    }
                  }
                }
                // For each OR gate.
                if ((!init_ID) && (!update_ID) && kmer_found &&
                    (exact_match || (!save_distances && !save_matches))) {
                  break;
                }
              }
              // For each k-mer.
              if ((!init_ID) && (!update_ID) && (!save_matches) && (!save_distances) && (num_matched >= c)) {
                break;
              }
              if (min_dist <= maximum_distance) {
                if (save_distances)
                  min_distances[min_dist]++;
                if (save_matches) {
                  match_cIDs.push_back(kmer_ID);
                  match_distances.push_back(min_dist);
                }
              }
            }
          }
          // For each cont.
          if ((!init_ID) && (!update_ID) && (!save_matches) && (num_matched >= c) && (!save_distances)) {
            break;
          }
        }

        // Try reverse complement.
        if ((num_matched < c) || save_distances || save_matches) {
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
                string kmer_str;
                if (i == 0) {
                  kmer_str = token.substr(i, int(k));
                  const char *ckmer = kmer_str.c_str();
                  encode_kmer_reverse(ckmer, b_enc, b_sig);
                  i = k - 1;
                } else {
                  kmer_str = token.substr(i, 1);
                  const char *ckmer = kmer_str.c_str();
                  update_kmer_reverse(ckmer, b_enc, b_sig);
                }

                bool kmer_found = false;
                bool closest_match = false;
                bool exact_match = false;
                uint8_t min_dist = KMER_LENGTH;
                uint16_t kmer_ID;

                for (uint64_t funci = 0; funci < l; funci++) {
                  kmer_sig = encode_kmer_bits(b_sig, shifts[funci], grab_bits[funci]);

                  tag = (kmer_sig >> ((2 * h) - t)) & tag_mask;
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
                      encid = get_encid(b * f_tmp + b * s_tmp + enc, encid_arr);

                      uint32_t encoding_ix = sigs_arr[f_tmp * b + s_tmp * b + enc];
                      if (encid == 0) {
                        test_enc = encode_arr_0[encoding_ix];
                      } else {
                        test_enc = encode_arr_1[encoding_ix];
                      }

                      uint8_t dist = hd(b_enc, test_enc);
                      check_distance(dist, p, min_dist, kmer_found, closest_match, exact_match,
                                     num_matched_reverse);

                      if (update_ID && exact_match) {
#pragma omp critical
                        {
                          update_kmer_cID(cID_arr_0, cID_arr_1, count_arr_0, count_arr_1, seen_0, seen_1, encid,
                                          encoding_ix, filename_cID, taxonomy_lookup);
                        }
                      }
                      if (init_ID && exact_match) {
#pragma omp critical
                        {
                          update_kmer_count(count_arr_0, count_arr_1, seen_0, seen_1, encid, encoding_ix);
                        }
                      }
                      if (closest_match && save_matches) {
                        if (encid == 0) {
                          kmer_ID = cID_arr_0[encoding_ix];
                        } else {
                          kmer_ID = cID_arr_1[encoding_ix];
                        }
                      }
                      // For each signature pointed row.
                      if ((!init_ID) && (!update_ID) && kmer_found &&
                          (exact_match || (!save_distances && !save_matches))) {
                        break;
                      }
                    }
                  }
                  // For each OR gate.
                  if ((!init_ID) && (!update_ID) && kmer_found &&
                      (exact_match || (!save_distances && !save_matches))) {
                    break;
                  }
                }
                // For each k-mer.
                if ((!init_ID) && (!update_ID) && (!save_matches) && (!save_distances) && (num_matched >= c)) {
                  break;
                }
                if (min_dist <= maximum_distance) {
                  if (save_distances)
                    reverse_min_distances[min_dist]++;
                  if (save_matches) {
                    rc_match_taxIDs.push_back(kmer_ID);
                    rc_match_distances.push_back(min_dist);
                  }
                }
              }
            }
            // For each cont.
            if ((!init_ID) && (!update_ID) && (!save_matches) && (num_matched >= c) && (!save_distances)) {
              break;
            }
          }
        }

        if (save_distances) {
#ifdef READ_PARALLELISM
#pragma omp critical
#endif
          { output_distances(ofs_kmer_distances, name, min_distances, reverse_min_distances); }
        }

        if (save_matches) {
#ifdef READ_PARALLELISM
#pragma omp critical
#endif
          {
            output_matches(ofs_match_information, name, match_cIDs, match_distances, rc_match_taxIDs,
                           rc_match_distances, cID_to_taxID);
          }
        }

        if ((num_matched < c) && (num_matched_reverse < c)) {
          if (unclassified_out) {
#ifdef READ_PARALLELISM
#pragma omp critical
#endif
            { output_reads(ofs_reads_unclassified, name, orig_read, line_third, line_fourth); }
          }
        } else if ((num_matched >= c) || (num_matched_reverse >= c)) {
#ifdef READ_PARALLELISM
#pragma omp atomic
#endif
          reads_matched += 1;
          if (classified_out) {
#ifdef READ_PARALLELISM
#pragma omp critical
#endif
            { output_reads(ofs_reads_classified, name, orig_read, line_third, line_fourth); }
          }
        }
      }
    }

    cout << query_file_path << " " << num_lines_read << " " << reads_matched << endl;

    ifs_reads_query.close();

    if (unclassified_out)
      ofs_reads_unclassified.close();
    if (classified_out)
      ofs_reads_classified.close();
    if (save_distances)
      ofs_kmer_distances.close();
    if (save_matches)
      ofs_match_information.close();
  }

  if (update_ID || init_ID) {
    // Write cIDs-0.
    for (int m = 0; m < encf_chunks; m++) {
      FILE *wfcID;
      wfcID = fopen(str_map_cID[m].c_str(), "wb");
      if (!wfcID) {
        cout << "Cannot open file for k-mer cIDs in the library directory!" << endl;
        exit(1);
      }
      fwrite(cID_arr_0 + enc_chunk_cumcounts_0[m], sizeof(uint16_t), enc_chunk_counts_0[m], wfcID);
      fclose(wfcID);
    }
    // Write cIDs-1.
    for (int m = 0; m < encf_chunks; m++) {
      FILE *wfcID;
      wfcID = fopen(str_map_cIDe[m].c_str(), "wb");
      if (!wfcID) {
        cout << "Cannot open file for k-mer cIDs in the library directory!" << endl;
        exit(1);
      }
      fwrite(cID_arr_1 + enc_chunk_cumcounts_1[m], sizeof(uint16_t), enc_chunk_counts_1[m], wfcID);
      fclose(wfcID);
    }
  }
  if (init_ID) {
    // Write counts-0.
    for (int m = 0; m < encf_chunks; m++) {
      FILE *wfcount;
      wfcount = fopen(str_map_count[m].c_str(), "wb");
      if (!wfcount) {
        cout << "Cannot open file for count-0 in the library directory!" << endl;
        exit(1);
      }
      fwrite(count_arr_0 + enc_chunk_cumcounts_0[m], sizeof(uint16_t), enc_chunk_counts_0[m], wfcount);
      fclose(wfcount);
    }
    // Write counts-1.
    for (int m = 0; m < encf_chunks; m++) {
      FILE *wfcount;
      wfcount = fopen(str_map_counte[m].c_str(), "wb");
      if (!wfcount) {
        cout << "Cannot open file for count-1 in the library directory!" << endl;
        exit(1);
      }
      fwrite(count_arr_1 + enc_chunk_cumcounts_1[m], sizeof(uint16_t), enc_chunk_counts_1[m], wfcount);
      fclose(wfcount);
    }
  }

  // Remember to delete array when it is done.
  delete[] sigs_arr;
  delete[] tag_arr;
  delete[] encode_arr_0;
  delete[] cID_arr_0;
  delete[] encode_arr_1;
  delete[] cID_arr_1;
  delete[] encid_arr;

  end = chrono::steady_clock::now();
  cout << "Done matching for all. Time so far: " << chrono::duration_cast<chrono::seconds>(end - start).count()
       << " seconds." << endl;

  return 0;
}

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
    if ((strcmp(entry->d_name, "..") != 0) && (strcmp(entry->d_name, ".") != 0)) {
      userString.push_back(string(path) + "/" + entry->d_name);
    }
  }

  closedir(dir);
  return (userString);
}

uint64_t encode_kmer_bits(uint64_t val, vector<int8_t> shifts, vector<int8_t> bits_to_grab) {
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

uint8_t get_encid(uint64_t sind, uint8_t encid_arr[]) {
  uint64_t eind = sind >> 3;
  uint64_t ebit = sind % 8;

  uint8_t encid;
  encid = (encid_arr[eind] >> (7 - ebit)) & 1;

  return (encid);
}

void check_distance(uint8_t &dist, uint64_t &p, uint8_t &min_dist, bool &kmer_found, bool &closest_match,
                    bool &exact_match, uint8_t &num_matched) {
  if ((!kmer_found) && (dist <= p)) {
    kmer_found = true;
    num_matched += 1;
  }
  if (dist < min_dist) {
    min_dist = dist;
    closest_match = true;
  } else {
    closest_match = false;
  }
  exact_match = dist == 0;
}

void output_reads(ofstream &ofs_reads, string name, string orig_read, string line_third, string line_fourth) {
  ofs_reads << name << endl;
  ofs_reads << orig_read << endl;
  // Output separator and quality.
  ofs_reads << line_third << endl;
  ofs_reads << line_fourth << endl;
}

void output_distances(ofstream &ofs_kmer_distances, string name, map<uint64_t, uint64_t> min_distances,
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

void output_matches(ofstream &ofs_match_information, string name, vector<uint16_t> match_cIDs,
                    vector<uint16_t> match_distances, vector<uint16_t> rc_match_taxIDs,
                    vector<uint16_t> rc_match_distances, map<uint16_t, uint64_t> &cID_to_taxID) {
  ofs_match_information << name << endl;

  ofs_match_information << "--";
  for (int i = 0; i < match_cIDs.size(); ++i) {
    ofs_match_information << " " << cID_to_taxID[match_cIDs[i]] << ":" << match_distances[i];
  }
  ofs_match_information << endl;

  ofs_match_information << "rc";
  for (int i = 0; i < rc_match_taxIDs.size(); ++i) {
    ofs_match_information << " " << cID_to_taxID[rc_match_taxIDs[i]] << ":" << rc_match_distances[i];
  }
  ofs_match_information << endl;
}

map<uint64_t, uint64_t> init_distance_map(uint64_t maximum_distance) {
  map<uint64_t, uint64_t> distances;
  for (uint64_t i = 0; i <= maximum_distance; ++i) {
    distances[i] = 0;
  }
  return distances;
}

void read_filename_map(string filename, map<string, uint64_t> &filename_to_taxID) {
  ifstream fmap;
  fmap.open(filename);
  if (!fmap) {
    cout << "Cannot open file for filename to taxonomic ID map." << endl;
    exit(1);
  }
  string key;
  string value;
  while (fmap >> key >> value) {
    filename_to_taxID[key] = stoi(value);
  }
  fmap.close();
}

uint16_t getLCA(uint16_t cID_0, uint16_t cID_1, unordered_map<uint16_t, vector<uint16_t>> &taxonomy_lookup) {
  vector<uint16_t> p_curr = taxonomy_lookup[cID_0];
  vector<uint16_t> p_new = taxonomy_lookup[cID_1];
  uint16_t cID = 1;
  for (int i = 0; i < p_curr.size(); ++i) {
    if ((p_new[i + p_new.size() - p_curr.size()] != 0) && (p_curr[i] != 0)) {
      if (p_curr[i] == p_new[i + p_new.size() - p_curr.size()]) {
        cID = p_curr[i];
        break;
      }
    } else {
      if (p_new[i + p_new.size() - p_curr.size()] != 0) {
        cID = p_new[i + p_new.size() - p_curr.size()];
        break;
      } else if (p_curr[i] != 0) {
        cID = p_curr[i];
        break;
      }
    }
  }
  return cID;
}

void update_kmer_cID(uint16_t cID_arr_0[], uint16_t cID_arr_1[], uint16_t count_arr_0[],
                     uint16_t count_arr_1[], vector<bool> &seen_0, vector<bool> &seen_1, uint8_t encid,
                     uint32_t encoding_ix, uint16_t filename_cID,
                     unordered_map<uint16_t, vector<uint16_t>> &taxonomy_lookup) {
  float p_update;
  float w = 5.0;
  float s = 4.0;
  random_device device;
  mt19937 gen(device());
  if (encid == 0) {
    if (!seen_0[encoding_ix]) {
      seen_0[encoding_ix] = true;
      p_update = min(1.0, pow(1.0 / w, 2) + (s / max(s, s + (float)count_arr_0[encoding_ix] - w)));
      bernoulli_distribution btrial(p_update);
      bool update = btrial(gen);
      if (update) {
        if (cID_arr_0[encoding_ix] == 0) {
          cID_arr_0[encoding_ix] = filename_cID;
        } else if (cID_arr_0[encoding_ix] == 1) {
        } else {
          cID_arr_0[encoding_ix] = getLCA(cID_arr_0[encoding_ix], filename_cID, taxonomy_lookup);
        }
      }
    }
  } else {
    if (!seen_1[encoding_ix]) {
      seen_1[encoding_ix] = true;
      p_update = min(1.0, pow(1.0 / w, 2) + (s / max(s, s + (float)count_arr_1[encoding_ix] - w)));
      bernoulli_distribution btrial(p_update);
      bool update = btrial(gen);
      if (update) {
        if (cID_arr_1[encoding_ix] == 0) {
          cID_arr_1[encoding_ix] = filename_cID;
        } else if (cID_arr_1[encoding_ix] == 1) {
        } else {
          cID_arr_1[encoding_ix] = getLCA(cID_arr_1[encoding_ix], filename_cID, taxonomy_lookup);
        }
      }
    }
  }
}

void update_kmer_count(uint16_t count_arr_0[], uint16_t count_arr_1[], vector<bool> &seen_0,
                       vector<bool> &seen_1, uint8_t encid, uint32_t encoding_ix) {
  if (encid == 0) {
    if (!seen_0[encoding_ix]) {
      count_arr_0[encoding_ix]++;
      seen_0[encoding_ix] = true;
    }
  } else {
    if (!seen_1[encoding_ix]) {
      count_arr_1[encoding_ix]++;
      seen_1[encoding_ix] = true;
    }
  }
}

void read_taxonomy_lookup(string filepath, unordered_map<uint16_t, vector<uint16_t>> &taxonomy_lookup,
                          map<uint64_t, uint16_t> &taxID_to_cID, map<uint16_t, uint64_t> &cID_to_taxID) {
  ifstream ftable;
  ftable.open(filepath);
  if (!ftable) {
    cout << "Cannot open file for taxonomic lookup table." << endl;
    exit(1);
  }
  vector<uint64_t> taxID_vec;
  unordered_map<uint64_t, vector<uint64_t>> taxID_lookup;

  for (string line; getline(ftable, line);) {
    istringstream iss(line);
    string taxIDstr;
    getline(iss, taxIDstr, ' ');
    uint64_t taxID = stoi(taxIDstr);
    taxID_vec.push_back(taxID);

    vector<uint64_t> lineage;
    string next_taxIDstr;

    while (getline(iss, next_taxIDstr, ',')) {
      uint64_t next_taxID = stoi(next_taxIDstr);
      lineage.push_back(next_taxID);
    }
    taxID_lookup.insert({taxID, lineage});
  }
  if (taxID_vec.size() > TAX_ID_LIMIT) {
    cout << "Total number of taxonomic IDs encountered exceeds limit." << endl;
    exit(1);
  }

  taxID_to_cID[1] = 1;
  cID_to_taxID[1] = 1;
  taxID_to_cID[0] = 0;
  cID_to_taxID[0] = 0;

  sort(taxID_vec.begin(), taxID_vec.end());

  for (uint16_t i = 0; i < taxID_vec.size(); ++i) {
    taxID_to_cID[taxID_vec[i]] = i + 2;
    cID_to_taxID[i + 2] = taxID_vec[i];
  }

  for (const auto &kv : taxID_lookup) {
    uint16_t cID = taxID_to_cID[kv.first];
    vector<uint16_t> lineage;
    for (const auto &taxID : kv.second) {
      lineage.push_back(taxID_to_cID[taxID]);
    }
    taxonomy_lookup.insert({cID, lineage});
  }
}
