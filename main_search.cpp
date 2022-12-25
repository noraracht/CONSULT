#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <new>
#include <omp.h>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <string>
#include <sys/types.h>
#include <time.h>
#include <unordered_map>
#include <vector>

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define VERSION 18.0
#define SL 32

using namespace std;

// Prototypes.
vector<string> list_dir(const char *path);
uint8_t hd(uint64_t x, uint64_t y);
uint8_t get_encid(uint64_t sind, uint8_t enc_arr_id[]);
uint64_t encodekmer_bits_rev(const char *s, vector<int> pos);
uint64_t encodekmer_bits(uint64_t val, vector<int8_t> shifts,
                         vector<int8_t> bits_to_grab);

void encodekmer(const char s[], uint64_t &b, uint64_t &b_sig);
void encodekmer_rev(const char s[], uint64_t &b, uint64_t &b_sig);

void update_kmer(const char *s, uint64_t &b, uint64_t &b_sig);
void update_kmer_rev(const char *s, uint64_t &b, uint64_t &b_sig);

int main(int argc, char *argv[]) {
  auto start = chrono::steady_clock::now();

  // Display version number.
  cout << "v." << std::fixed << std::setprecision(1) << VERSION << endl;

  char *input_library_dir = NULL;
  char *query_dir = NULL;
  uint64_t c_value = 1;
  uint64_t thread_count = 1;
  bool save_distances = false;

  int cf_tmp;
  opterr = 0;

  while (1) {
    static struct option long_options[] = {
        {"input-library-dir", 1, 0, 'i'}, {"query-dir", 1, 0, 'q'},
        {"c-value", 1, 0, 'c'},           {"thread-count", 1, 0, 't'},
        {"save-distances", 0, 0, 's'},    {0, 0, 0, 0},
    };

    int option_index = 0;
    cf_tmp = getopt_long(argc, argv, "i:q:c:t:s", long_options, &option_index);

    if ((optarg != NULL) && (*optarg == '-')) {
      cf_tmp = ':';
    }

    if (cf_tmp == -1)
      break;
    switch (cf_tmp) {
    case 'i':
      input_library_dir = optarg;
      break;
    case 'c':
      c_value = max(atoi(optarg), 1); // Default is 1.
      break;
    case 't':
      thread_count = max(atoi(optarg), 1); // Default is 1.
      break;
    case 'q':
      query_dir = optarg;
      break;
    case 's':
      save_distances = true;
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
        fprintf(stderr, "Unknown option `-%c'.\n", optopt);
      else
        fprintf(stderr, "Unknown option '%s'.\n", argv[optind - 1]);
      return 1;
    default:
      abort();
    }
  }

  if (argc <= 10) {
    printf("-- Arguments supplied are \nMap %s\nQuery %s\nThreads %ld\n",
           input_library_dir, query_dir, thread_count);
  } else if (argc > 10) {
    printf("Too many arguments are supplied.\n");
    exit(0);
  }

  // Read map.
  string map_name = input_library_dir;
  string map_meta = map_name + "_meta";
  string path = map_name + "/" + map_meta;

  FILE *fmeta = fopen(path.c_str(), "rb");
  if (!fmeta) {
    cout << "Cannot open file!" << endl;
    return 1;
  }

  // Set confidence threshold, should be at 0.
  uint64_t c = c_value;
  // Read parameters from input file.
  uint64_t p;
  uint64_t L;
  float alpha;
  uint64_t K;
  uint64_t sigs_arr_size;
  uint64_t new_tag_arr_size;
  uint64_t kmer_count;
  uint64_t encli_0;
  uint64_t encli_1;
  uint64_t enc_arr_id_size;

  fread(&p, sizeof(uint64_t), 1, fmeta);
  fread(&L, sizeof(uint64_t), 1, fmeta);
  fread(&alpha, sizeof(float), 1, fmeta);
  fread(&K, sizeof(uint64_t), 1, fmeta);
  fread(&sigs_arr_size, sizeof(uint64_t), 1, fmeta);
  fread(&new_tag_arr_size, sizeof(uint64_t), 1, fmeta);
  fread(&kmer_count, sizeof(uint64_t), 1, fmeta);
  fread(&encli_0, sizeof(uint64_t), 1, fmeta);
  fread(&encli_1, sizeof(uint64_t), 1, fmeta);
  fread(&enc_arr_id_size, sizeof(uint64_t), 1, fmeta);

  string line;
  uint64_t l;

  cout << "c = " << c << endl;
  cout << "p = " << p << '\n';
  cout << "SL = " << SL << '\n';
  cout << "L = " << L << '\n';
  cout << "alpha = " << alpha << '\n';
  cout << "Using K = " << K << '\n';
  cout << "k-mer count = " << kmer_count << endl;
  cout << "k-mer count array 0  = " << encli_0 << endl;
  cout << "k-mer count array 1  = " << encli_1 << endl;
  cout << "Tag array size  = " << new_tag_arr_size << endl;
  cout << "Signature array size = " << sigs_arr_size << endl;
  cout << "Encoding id array size  = " << enc_arr_id_size << endl;

  int vec_size = 0;
  vector<vector<int8_t>> shifts;
  vector<vector<int8_t>> grab_bits;

  for (l = 0; l < L; l++) {
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
  uint64_t SIGS_COL_COUNT;
  uint64_t SIGS_ROW_COUNT;
  fread(&SIGS_COL_COUNT, sizeof(uint64_t), 1, fmeta);
  fread(&SIGS_ROW_COUNT, sizeof(uint64_t), 1, fmeta);

  // Read tag size, partition count, tag mask and big sig mask.
  uint64_t tag_size;
  fread(&tag_size, sizeof(uint64_t), 1, fmeta);
  uint64_t partitions;
  fread(&partitions, sizeof(uint64_t), 1, fmeta);
  uint64_t tag_mask;
  fread(&tag_mask, sizeof(uint64_t), 1, fmeta);
  uint64_t big_sig_mask;
  fread(&big_sig_mask, sizeof(uint64_t), 1, fmeta);

  // Display map information.
  cout << "Sigs row count = " << SIGS_ROW_COUNT << '\n';
  cout << "Columns = " << unsigned(SIGS_COL_COUNT) << endl;
  cout << "Partitions = " << unsigned(partitions) << endl;
  cout << "Tag size = " << unsigned(tag_size) << '\n';
  cout << "Tag mask = " << unsigned(tag_mask) << '\n';
  cout << "Big sig mask = " << big_sig_mask << '\n';

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
    cout << "-- Done sigs allocation" << endl;
  } catch (std::bad_alloc &ba) {
    std::cerr << "Exception 'bad_alloc' is caught: " << ba.what() << endl;
  }

  // Allocate tag array.
  int8_t *tag_arr;

  try {
    tag_arr = new int8_t[new_tag_arr_size];
    cout << "-- Done tag allocation" << endl;
  } catch (std::bad_alloc &ba) {
    std::cerr << "Exception bad_alloc is caught: " << ba.what() << endl;
  }

  // Allocate and read encoding array.
  uint64_t *encode_arr_0;
  uint64_t *encode_arr_1;

  try {
    encode_arr_0 = new uint64_t[encli_0];
    cout << "-- Done encoding array 0 allocation" << endl;
  } catch (std::bad_alloc &ba) {
    std::cerr << "Exception bad_alloc is caught: " << ba.what() << endl;
  }

  try {
    encode_arr_1 = new uint64_t[encli_1];
    cout << "-- Done encoding array 1 allocation" << endl;
  } catch (std::bad_alloc &ba) {
    std::cerr << "Exception bad_alloc is caught: " << ba.what() << endl;
  }

  // Allocate enc array id array.
  uint8_t *enc_arr_id;

  try {
    enc_arr_id = new uint8_t[enc_arr_id_size];
    cout << "-- Done enc id allocation" << endl;
  } catch (std::bad_alloc &ba) {
    std::cerr << "Exception bad_alloc is caught: " << ba.what() << '\n';
  }

  uint64_t total_sigs_read = 0;
  uint64_t total_tags_read = 0;
  uint64_t num_pairs_0 = 0;
  uint64_t num_pairs_1 = 0;
  uint64_t total_encid_read = 0;

  // Read files in parallel.
  vector<string> str_map_sig;
  vector<string> str_map_tag;
  vector<string> str_map_enc;
  vector<string> str_map_ence;
  vector<string> str_map_encid;

  for (int m = 0; m < sigf_chunks; m++) {
    string map_sig = map_name + "_sig" + to_string(m);
    path = map_name + "/" + map_sig;
    str_map_sig.push_back(path);
  }
  for (int m = 0; m < tagf_chunks; m++) {
    string map_tag = map_name + "_tag" + to_string(m);
    path = map_name + "/" + map_tag;
    str_map_tag.push_back(path);
  }
  for (int m = 0; m < encf_chunks; m++) {
    string map_enc = map_name + "_enc" + to_string(m);
    path = map_name + "/" + map_enc;
    str_map_enc.push_back(path);
  }
  for (int m = 0; m < encf_chunks; m++) {
    string map_ence = map_name + "_ence" + to_string(m);
    path = map_name + "/" + map_ence;
    str_map_ence.push_back(path);
  }
  for (int m = 0; m < sigf_chunks; m++) {
    string map_encid = map_name + "_encid" + to_string(m);
    path = map_name + "/" + map_encid;
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
        cout << "Cannot open file!" << endl;
        exit(0);
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
        cout << "Cannot open file!" << endl;
        exit(0);
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
        cout << "Cannot open file!" << endl;
        exit(0);
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
        cout << "Cannot open file!" << endl;
        exit(0);
      }

      // Read encodings.
      temp_count = fread(encode_arr_1 + enc_chunk_cumcounts_1[m],
                         sizeof(uint64_t), enc_chunk_counts_1[m], fence);
#pragma omp atomic

      num_pairs_1 += temp_count;

      fclose(fence);
    }
  }

// Read encid array.
#pragma omp parallel num_threads(thread_count)                                 \
    shared(enc_arr_id, total_encid_read)
  {
#pragma omp for
    for (int m = 0; m < sigf_chunks; m++) {
      FILE *fencid;
      fencid = fopen(str_map_encid[m].c_str(), "rb");

      if (!fencid) {
        cout << "Cannot open file!" << endl;
        exit(0);
      }

      // Read sigs.
      temp_count = fread(enc_arr_id + encid_chunk_cumcounts[m], sizeof(uint8_t),
                         encid_chunk_counts[m], fencid);
#pragma omp atomic

      total_encid_read += temp_count;

      fclose(fencid);
    }
  }

  cout << "Sigs read:  " << total_sigs_read << endl;
  cout << "Tags read:  " << total_tags_read << endl;
  cout << "Encodings 0 read: " << num_pairs_0 << endl;
  cout << "Encodings 1 read:" << num_pairs_1 << endl;
  cout << "Enc id read: " << total_encid_read << endl;

  auto end = chrono::steady_clock::now();
  cout << "-- Done reading. Now matching. Time so far: "
       << chrono::duration_cast<chrono::seconds>(end - start).count()
       << " seconds" << endl;

  // Read input fastq.
  const char *dir = query_dir;

  vector<string> file_list;
  if (dir == NULL)
    return (1);
  else
    file_list = list_dir(dir);

  int file_count = file_list.size();

#pragma omp parallel num_threads(thread_count)
  {
#pragma omp for schedule(dynamic)
    for (uint64_t f = 0; f < file_count; f++) {

      string input_fq = file_list[f];

      stringstream stream;
      stream << std::fixed << std::setprecision(2) << alpha;
      std::string alpha_s = stream.str();

      string input_fq_truct = input_fq.substr(input_fq.find_last_of("/") + 1);
      string output_fname = "ucseq_" + input_fq_truct;

      ifstream ifs(input_fq);

      ofstream outputFile;
      outputFile.open(output_fname);

      uint64_t lines_read = 0;
      uint64_t reads_matched = 0;

      string line_of_file;
      string name;
      string token;

      uint64_t b;
      uint64_t b_sig;

      uint64_t test_enc;
      uint8_t enc_arr_ind;

      int8_t tag;
      uint64_t kmer_sig;
      uint64_t big_sig_hash;
      uint64_t enc_start;
      uint64_t enc_end;

      while (!ifs.eof()) {
        getline(ifs, line_of_file);

        if (!ifs)
          break;

        if (lines_read % 4 == 0) {
          name = line_of_file;
        } else if (lines_read % 4 == 1) {
          int matched = 0;

          string line_of_file_orig = line_of_file;

          istringstream iss(line_of_file);

          while (getline(iss, token, 'N')) {
            if (token.length() >= int(SL)) {
              b = 0;
              b_sig = 0;

              for (uint64_t i = 0; i < token.length(); i++) {
                if (i == 0) {
                  string kmer_str = token.substr(i, int(SL));
                  const char *ckmer = kmer_str.c_str();
                  encodekmer(ckmer, b, b_sig);
                  i = SL - 1;

                } else {
                  string kmer_str = token.substr(i, 1);
                  const char *ckmer = kmer_str.c_str();
                  update_kmer(ckmer, b, b_sig);
                }

                bool kmerfound = false;

                for (int64_t funci = 0; funci < L; funci++) {
                  kmer_sig =
                      encodekmer_bits(b_sig, shifts[funci], grab_bits[funci]);

                  // Get first 2 bits of signature (effectively bits 28 - 27) of
                  // 32 bit encoding as partition numbers.
                  tag = (kmer_sig >> ((2 * K) - tag_size)) & tag_mask;

                  // Get last 26 bits of signature (effectively bits 26 -1 ) of
                  // 32 bit encoding as sigs row number.
                  big_sig_hash = kmer_sig & big_sig_mask;

                  if (tag_arr[SIGS_ROW_COUNT * partitions * funci +
                              big_sig_hash * partitions + tag] != -1) {
                    if (tag == 0) {
                      enc_start = 0;
                    } else {
                      enc_start = tag_arr[SIGS_ROW_COUNT * partitions * funci +
                                          big_sig_hash * partitions + tag - 1];
                    }
                    enc_end = tag_arr[SIGS_ROW_COUNT * partitions * funci +
                                      big_sig_hash * partitions + tag];

                    for (uint64_t enc = enc_start; enc < enc_end; enc++) {
                      // Set encoding array id.
                      enc_arr_ind = get_encid(
                          SIGS_ROW_COUNT * SIGS_COL_COUNT * partitions * funci +
                              big_sig_hash * SIGS_COL_COUNT * partitions + enc,
                          enc_arr_id);

                      if (enc_arr_ind == 0) {
                        test_enc = encode_arr_0
                            [sigs_arr[SIGS_ROW_COUNT * SIGS_COL_COUNT *
                                          partitions * funci +
                                      big_sig_hash * SIGS_COL_COUNT *
                                          partitions +
                                      enc]];
                      } else {
                        test_enc = encode_arr_1
                            [sigs_arr[SIGS_ROW_COUNT * SIGS_COL_COUNT *
                                          partitions * funci +
                                      big_sig_hash * SIGS_COL_COUNT *
                                          partitions +
                                      enc]];
                      }

                      int8_t dist = hd(b, test_enc);
                      kmerfound = dist <= p;
                      // For each signature pointed row.
                      if (kmerfound) {
                        break;
                      }
                    }
                  }
                  // For each OR gate.
                  if (kmerfound) {
                    matched += 1;
                    break;
                  }
                }
                // For each k-mer.
                if (matched >= c) {
                  break;
                }
              }
            }
            // For each line.
            if (matched >= c) {
              break;
            }
          }

          // Try reverse complement.
          if (matched < c) {
            matched = 0;

            int len = strlen(line_of_file.c_str());
            char swap;

            for (int i = 0; i < len / 2; i++) {
              swap = line_of_file[i];
              line_of_file[i] = line_of_file[len - i - 1];
              line_of_file[len - i - 1] = swap;
            }

            istringstream iss(line_of_file);

            while (getline(iss, token, 'N')) {
              if (token.length() >= int(SL)) {

                b = 0;
                b_sig = 0;

                for (uint64_t i = 0; i < token.length(); i++) {
                  if (i == 0) {
                    string kmer_str = token.substr(i, int(SL));
                    const char *ckmer = kmer_str.c_str();
                    encodekmer_rev(ckmer, b, b_sig);
                    i = SL - 1;
                  } else {
                    string kmer_str = token.substr(i, 1);
                    const char *ckmer = kmer_str.c_str();
                    update_kmer_rev(ckmer, b, b_sig);
                  }

                  bool kmerfound = false;

                  for (uint64_t funci = 0; funci < L; funci++) {
                    kmer_sig =
                        encodekmer_bits(b_sig, shifts[funci], grab_bits[funci]);

                    // Get first 2 bits of signature (effectively bits 28 - 27)
                    // of 32 bit encoding as partition numbers.
                    tag = (kmer_sig >> ((2 * K) - tag_size)) & tag_mask;

                    // Get last 26 bits of signature (effectively bits 26 -1 )
                    // of 32 bit encoding as sigs row number.
                    big_sig_hash = kmer_sig & big_sig_mask;

                    if (tag_arr[SIGS_ROW_COUNT * partitions * funci +
                                big_sig_hash * partitions + tag] != -1) {
                      if (tag == 0) {
                        enc_start = 0;
                      } else {
                        enc_start =
                            tag_arr[SIGS_ROW_COUNT * partitions * funci +
                                    big_sig_hash * partitions + tag - 1];
                      }
                      enc_end = tag_arr[SIGS_ROW_COUNT * partitions * funci +
                                        big_sig_hash * partitions + tag];

                      for (uint64_t enc = enc_start; enc < enc_end; enc++) {
                        enc_arr_ind = get_encid(
                            SIGS_ROW_COUNT * SIGS_COL_COUNT * partitions *
                                    funci +
                                big_sig_hash * SIGS_COL_COUNT * partitions +
                                enc,
                            enc_arr_id);

                        if (enc_arr_ind == 0) {
                          test_enc = encode_arr_0
                              [sigs_arr[SIGS_ROW_COUNT * SIGS_COL_COUNT *
                                            partitions * funci +
                                        big_sig_hash * SIGS_COL_COUNT *
                                            partitions +
                                        enc]];
                        } else {
                          test_enc = encode_arr_1
                              [sigs_arr[SIGS_ROW_COUNT * SIGS_COL_COUNT *
                                            partitions * funci +
                                        big_sig_hash * SIGS_COL_COUNT *
                                            partitions +
                                        enc]];
                        }

                        int8_t dist = hd(b, test_enc);
                        kmerfound = dist <= p;
                        // For each signature pointed row.
                        if (kmerfound) {
                          break;
                        }
                      }
                    }
                    // For each OR gate.
                    if (kmerfound) {
                      matched += 1;
                      break;
                    }
                  }
                  // For each k-mer.
                  if (matched >= c) {
                    break;
                  }
                }
              }
              // For each line.
              if (matched >= c) {
                break;
              }
            }
          }

          if (matched < c) {
            outputFile << name << endl;
            outputFile << line_of_file_orig << endl;

            // Output separator and quality.
            getline(ifs, line_of_file);
            ++lines_read;
            outputFile << line_of_file << endl;

            getline(ifs, line_of_file);
            ++lines_read;
            outputFile << line_of_file << endl;

          } else if (matched >= c) {
            reads_matched += 1;
          }
        } else {
        }
        ++lines_read;
      }

#pragma omp critical
      { cout << input_fq << " " << lines_read << " " << reads_matched << endl; }

      ifs.close();
      outputFile.close();
    }
  }

  // Remember to delete array when it is done.
  delete[] sigs_arr;
  delete[] tag_arr;
  delete[] encode_arr_0;
  delete[] encode_arr_1;
  delete[] enc_arr_id;

  end = chrono::steady_clock::now();
  cout << "-- Done matching for all. Time so far: "
       << chrono::duration_cast<chrono::seconds>(end - start).count()
       << " seconds" << endl;

  return 0;
}

// Function definitions.

void encodekmer(const char *s, uint64_t &b, uint64_t &b_sig) {
  for (int i = 0; i < int(SL); i++) {
    b = b << 1;
    b_sig = b_sig << 2;

    if (s[i] == 'T') {
      b += 4294967297;
      b_sig += 3;
    } else if (s[i] == 'G') {
      b += 4294967296;
      b_sig += 2;
    } else if (s[i] == 'C') {
      b += 1;
      b_sig += 1;
    } else {
      b += 0;
      b_sig += 0;
    }
  }
}

void encodekmer_rev(const char *s, uint64_t &b, uint64_t &b_sig) {
  for (int i = 0; i < int(SL); i++) {
    b = b << 1;
    b_sig = b_sig << 2;

    if (s[i] == 'A') {
      b += 4294967297;
      b_sig += 3;
    } else if (s[i] == 'C') {
      b += 4294967296;
      b_sig += 2;
    } else if (s[i] == 'G') {
      b += 1;
      b_sig += 1;
    } else {
      b += 0;
      b_sig += 0;
    }
  }
}

void update_kmer(const char *s, uint64_t &b, uint64_t &b_sig) {
  b = b << 1;
  b_sig = b_sig << 2;

  // Create a mask that has set 32 bit only and set bit 32 to 0.
  uint64_t mask = 4294967297;
  b = b & ~mask;

  if (s[0] == 'T') {
    b += 4294967297;
    b_sig += 3;
  } else if (s[0] == 'G') {
    b += 4294967296;
    b_sig += 2;
  } else if (s[0] == 'C') {
    b += 1;
    b_sig += 1;
  } else {
    b += 0;
    b_sig += 0;
  }
}

void update_kmer_rev(const char *s, uint64_t &b, uint64_t &b_sig) {
  b = b << 1;
  b_sig = b_sig << 2;

  // Create a mask that has set 32 bit only and set bit 32 to 0.
  uint64_t mask = 4294967297;
  b = b & ~mask;

  if (s[0] == 'A') {
    b += 4294967297;
    b_sig += 3;
  } else if (s[0] == 'C') {
    b += 4294967296;
    b_sig += 2;
  } else if (s[0] == 'G') {
    b += 1;
    b_sig += 1;
  } else {
    b += 0;
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

uint64_t encodekmer_bits(uint64_t val, vector<int8_t> shifts,
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

uint64_t encodekmer_bits_rev(const char *s, vector<int> pos) {
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

uint8_t get_encid(uint64_t sind, uint8_t enc_arr_id[]) {
  uint64_t eind = sind >> 3;
  uint64_t ebit = sind % 8;

  uint8_t enc_arr_ind;
  enc_arr_ind = (enc_arr_id[eind] >> (7 - ebit)) & 1;

  return (enc_arr_ind);
}
