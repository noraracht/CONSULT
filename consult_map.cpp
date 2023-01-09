#include <algorithm>           // for sort, count, fill
#include <bits/getopt_core.h>  // for optarg, optopt, optind, opterr
#include <bits/stdint-uintn.h> // for uint64_t, uint8_t, uint32_t, uint16_t
#include <chrono>              // for seconds, duration_cast, operator-
#include <cmath>               // for round, pow, ceil
#include <cstdint>             // for int8_t, uint64_t, uint8_t
#include <cstdio>              // for fwrite, fclose, fopen, fprintf, FILE
#include <cstdlib>
#include <cstring> // for strerror
#include <ctype.h> // for isprint
#include <ctype.h>
#include <errno.h>            // for errno
#include <ext/alloc_traits.h> // for __alloc_traits<>::value_type
#include <fstream>
#include <functional> // for greater
#include <getopt.h>   // for getopt_long, option
#include <iomanip>    // for operator<<, setprecision
#include <iostream>   // for operator<<, endl, basic_ostream, ostream
#include <memory>     // for allocator_traits<>::value_type
#include <new>        // for bad_alloc
#include <numeric>    // for accumulate
#include <sstream>
#include <stdio.h>
#include <stdlib.h>   // for atoi, exit, rand, srand, abort
#include <string>     // for operator+, string, allocator, char_tr...
#include <sys/stat.h> // for mkdir, S_IROTH, S_IRWXG, S_IRWXU, S_I...
#include <sys/types.h>
#include <time.h> // for time
#include <tuple>  // for operator<, get, swap, make_tuple, tuple
#include <unistd.h>
#include <vector> // for vector, vector<>::value_type

#define VERSION 17.1
#define KMER_LENGTH 32

#define SIGF_CHUNKS 24
#define TAGF_CHUNKS 24
#define ENCF_CHUNKS 24

using namespace std;

// Prototypes.
uint64_t encodekmer_bits(uint64_t val, vector<int8_t> shifts,
                         vector<int8_t> bits_to_grab);

void encodekmer(const char s[], uint64_t &b, uint64_t &b_sig);

uint64_t file_read(istream &is, vector<char> &buff);
uint64_t count_lines(const vector<char> &buff, int size);

int main(int argc, char *argv[]) {
  auto start = chrono::steady_clock::now();
  srand(time(NULL));

  // Display version number.
  cout << "v." << std::fixed << std::setprecision(1) << VERSION << endl;

  string input_fasta_file = string();
  string output_library_dir = string();

  uint16_t k = KMER_LENGTH;
  /* Defaults */
  uint64_t p = 3;
  uint64_t l = 2;
  uint64_t h = 15;
  /* float alpha = 0.95; // This is just to determine correct h value. */
  /* uint64_t h = round(log(1 - (pow((1 - alpha), (1 / float(l))))) / */
  /*                    log(1 - float(p) / float(k))); */

  uint8_t t = 2;                                 // Tag size in bits.
  uint64_t partitions = pow(2, t);               // # of partitions; 2^(t)
  uint64_t sigs_col_count = 7;                   // # of columns per, i.e., b.
  uint64_t sigs_row_count = pow(2, (2 * h) - t); // # of rows; 2^(2h-t)

  int cf_tmp;
  opterr = 0;

  while (1) {
    static struct option long_options[] = {
        {"input-fasta-file", 1, 0, 'i'},
        {"output-library-directory", 1, 0, 'o'},
        {"p-value", 1, 0, 'p'},
        {"l-value", 1, 0, 'l'},
        {"h-value", 1, 0, 'h'},
        {"tag-size", 1, 0, 't'},
        {"column-per-tag", 1, 0, 'b'},
        {0, 0, 0, 0},
    };

    int option_index = 0;
    cf_tmp =
        getopt_long(argc, argv, "i:o:p:l:h:t:b:", long_options, &option_index);

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
      output_library_dir = optarg;
      break;
    case 'p':
      p = atoi(optarg); // Default value is 3.
      break;
    case 'l':
      l = atoi(optarg); // Default value is 2.
      break;
    case 'h':
      h = atoi(optarg); // Default value is 15.
      break;
    case 't':
      // Number of partitions is 2^t.
      t = atoi(optarg); // Default value is 2.
      break;
    case 'b':
      // Total number of columns is partitions * b.
      sigs_col_count = atoi(optarg); // Default value is 7.
      break;
    case ':':
      printf("Missing option for '-%s'.\n", argv[optind - 2]);
      if (long_options[option_index].has_arg == 1) {
        return 1;
      }
      break;
    case '?':
      if (optopt == 'i')
        fprintf(stderr, "Option '-%c' requires an argument.\n", optopt);
      else if (optopt == 'o')
        fprintf(stderr, "Option '-%c' requires an argument.\n", optopt);
      else if (isprint(optopt))
        fprintf(stderr, "Unknown option '-%c'.\n", optopt);
      else
        fprintf(stderr, "Unknown option '%s'.\n", argv[optind - 1]);
      return 1;
    default:
      abort();
    }
  }

  uint64_t kmer_count = 0;

  if (argc <= 15) {
    cout << "-- Arguments supplied are;" << endl;
    cout << "Input file : " << input_fasta_file << endl;
    cout << "Library directory : " << output_library_dir << endl;

    const int SZ = 1024 * 1024;
    vector<char> buff(SZ);

    // Open input file.
    ifstream ifs(input_fasta_file);
    if (!ifs.is_open()) {
      std::cout << "Cannot open file!" << endl;
      exit(0);
    }
    while (uint64_t cc = file_read(ifs, buff)) {
      kmer_count += count_lines(buff, cc);
    }
    ifs.close();
  } else {
    printf("Too many arguments are supplied.\n");
    exit(0);
  }

  ifstream ifs(input_fasta_file);
  if (!ifs.is_open()) {
    std::cout << "Cannot open file!" << endl;
    exit(0);
  }

  cout << "k-mer count = " << kmer_count << endl;
  cout << "k = " << k << '\n';
  cout << "p = " << p << '\n';
  cout << "l = " << l << '\n';
  cout << "Using h = " << h << '\n';
  /* cout << "alpha = " << alpha << '\n'; */

  // Allocate signature array.
  uint32_t *sigs_arr;
  uint64_t sigs_arr_size = sigs_row_count * sigs_col_count * partitions * l;

  try {
    sigs_arr = new uint32_t[sigs_arr_size];
    cout << "-- Done sigs allocation" << endl;
  } catch (std::bad_alloc &ba) {
    std::cerr << "Exception bad_alloc is caught: " << ba.what() << endl;
  }

  // Allocate sigs indicator array.
  uint8_t *sigs_indicator_arr;
  uint64_t sigs_indicator_arr_size = sigs_row_count * l;

  try {
    sigs_indicator_arr = new uint8_t[sigs_indicator_arr_size];
    cout << "-- Done indicator allocation" << endl;
  } catch (std::bad_alloc &ba) {
    std::cerr << "Exception bad_alloc is caught: " << ba.what() << endl;
  }

  // Initialize indicator array to 0.
  for (uint64_t enc = 0; enc < sigs_indicator_arr_size; enc++) {
    sigs_indicator_arr[enc] = 0;
  }

  // Allocate tag array.
  int8_t *tag_arr;
  uint64_t tag_arr_size = sigs_row_count * sigs_col_count * partitions * l;

  try {
    tag_arr = new int8_t[tag_arr_size];
    cout << "-- Done tag allocation" << endl;
  } catch (std::bad_alloc &ba) {
    std::cerr << "Exception bad_alloc is caught: " << ba.what() << endl;
  }

  /* Initialize tag array to -1. */
  // for (uint64_t j = 0; j < sigs_row_count * sigs_col_count * partitions * l;
  //      j++) {
  //   tag_arr[j] = -1;
  // }

  uint64_t MAX_ENC_CNT = 4000000000;

  // Allocate encoding array.
  uint64_t *encode_arr_0; // AT
  uint64_t *encode_arr_1; // CG
  uint64_t encode_arr_size;

  if (kmer_count < MAX_ENC_CNT) {
    encode_arr_size = kmer_count;
  } else {
    encode_arr_size = MAX_ENC_CNT + 5;
  }

  try {
    encode_arr_0 = new uint64_t[encode_arr_size];
    cout << "-- Done encoding array 0 allocation " << endl;
  } catch (std::bad_alloc &ba) {
    std::cerr << "Exception bad_alloc is caught: " << ba.what() << endl;
  }

  try {
    encode_arr_1 = new uint64_t[encode_arr_size];
    cout << "-- Done encoding array 1 allocation " << endl;
  } catch (std::bad_alloc &ba) {
    std::cerr << "Exception bad_alloc is caught: " << ba.what() << endl;
  }

  // Allocate enc array id array.
  uint8_t *enc_arr_id;
  uint64_t enc_arr_id_size = sigs_arr_size;

  try {
    enc_arr_id = new uint8_t[enc_arr_id_size];
    cout << "-- Done enc id array allocation " << endl;
  } catch (std::bad_alloc &ba) {
    std::cerr << "Exception bad_alloc is caught: " << ba.what() << endl;
  }

  // Initialize enc_arr_id to 0.
  for (uint64_t j = 0; j < enc_arr_id_size; j++) {
    enc_arr_id[j] = 0;
  }

  // Create a directory (Linux dependent functionality).
  const int dir_err =
      mkdir(output_library_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  if (-1 == dir_err) {
    cerr << "Cannot create directory! Error: " << strerror(errno) << endl;
    exit(1);
  }

  // Generate mask.
  uint64_t m, n;
  uint64_t included_kmers_counter[l];
  vector<vector<int>> positions(l, vector<int>(h));

  srand(time(NULL));

  for (int i = 0; i < l; i++) {
    vector<int> rand_num;

    // Initialize lost k-mers array.
    included_kmers_counter[i] = 0;

    for (m = 0; m < h; m++) {
      n = rand() % k;
      if (count(rand_num.begin(), rand_num.end(), n)) {
        m -= 1;
      } else {
        rand_num.push_back(n);
      }
    }

    cout << "Positions for l = " << i << " >> ";
    sort(rand_num.begin(), rand_num.end(), std::greater<int>());
    for (int j = 0; j < (h - 1); j++) {
      positions[i][j] = rand_num[j];
      cout << positions[i][j] << ", ";
    }
    cout << positions[i][h - 1] << endl;

    rand_num.clear();
  }

  vector<vector<int8_t>> shifts;
  vector<vector<int8_t>> grab_bits;

  for (int i = 0; i < l; i++) {
    // Construct a vector of integers.
    vector<int8_t> v;
    vector<int8_t> g;
    int lp = 31;
    int jp = 0;

    for (int j = 0; j < h; j++) {
      if (j == 0) {
        v.push_back((lp - positions[i][j]) * 2);
        lp = positions[i][j];
        jp += 2;

      } else if (positions[i][j - 1] - positions[i][j] != 1) {
        v.push_back((lp - positions[i][j]) * 2);
        lp = positions[i][j];
        g.push_back(jp);
        jp = 2;
      } else {
        jp += 2;
      }
    }
    g.push_back(jp);
    // Push back above one-dimensional vector.
    v.push_back(-1);
    g.push_back(-1);

    shifts.push_back(v);
    grab_bits.push_back(g);
  }

  string line;
  uint64_t li = 0;
  uint64_t encli_0 = 0;
  uint64_t encli_1 = 0;

  uint64_t b;
  uint64_t b_sig;

  bool kmer_written;

  int8_t tag;
  uint64_t sig_hash;
  uint64_t big_sig_hash;

  uint8_t enc_array_ind = 0;

  // Compute mask values.
  uint64_t tag_mask = 0;
  for (int i = 0; i < t; i++) {
    tag_mask += (pow(2, i)); // 111 if tag = 3
  }
  uint64_t big_sig_mask = 0;
  for (int i = 0; i < ((2 * h) - t); i++) {
    big_sig_mask += (pow(2, i));
  }

  while (std::getline(ifs, line)) {
    if (line.rfind('>', 0) != 0) {
      const char *cline = line.c_str();
      b = 0;
      b_sig = 0;
      kmer_written = false;

      encodekmer(cline, b, b_sig);

      if (li % 1000000 == 0) {
        cout << "-- Encoding " << li << endl;
      }

      // Set encoding array id.
      enc_array_ind = rand() % 2;

      // Check if max # of encodings (0xFFFFFFFF or 4294967295) isn't exceeded.
      if ((encli_0 >= MAX_ENC_CNT) || (encli_1 >= MAX_ENC_CNT)) {
        break;
      }

      for (int i = 0; i < l; i++) {
        sig_hash = encodekmer_bits(b_sig, shifts[i], grab_bits[i]);

        // Get first 2 bits of signature (effectively bits 28 - 27) of 32 bit
        // encoding as tag.
        tag = (sig_hash >> ((2 * h) - t)) & tag_mask;

        // Get last 26 bits of signature (effectively bits 26 - 1 ) of 32 bit
        // encoding as sigs row number.
        big_sig_hash = sig_hash & big_sig_mask;

        uint64_t tmp_idx =
            (sigs_row_count * sigs_col_count * partitions * i) +
            (big_sig_hash * sigs_col_count * partitions) +
            sigs_indicator_arr[sigs_row_count * i + big_sig_hash];

        // Check is row space is available for forward k-mer.
        if (sigs_indicator_arr[sigs_row_count * i + big_sig_hash] <
            sigs_col_count * partitions) {
          // Populate tag array.
          tag_arr[tmp_idx] = tag;

          if (enc_array_ind == 0) { // First letter is either A or C.
            // Populate sigs array.
            sigs_arr[tmp_idx] = encli_0;

          } else if (enc_array_ind == 1) // First letter is either G or T.
          {
            // Populate sigs array.
            sigs_arr[tmp_idx] = encli_1;
          }

          // Populate enc array id array.
          enc_arr_id[tmp_idx] = enc_array_ind;

          // Increment indicator array.
          sigs_indicator_arr[sigs_row_count * i + big_sig_hash] += 1;

          // Increment k-mer counter.
          included_kmers_counter[i] += 1;

          kmer_written = true;
        }
      }

      if (kmer_written) {
        if (enc_array_ind == 0) {
          encode_arr_0[encli_0] = b;
          encli_0 += 1;
        } else {
          encode_arr_1[encli_1] = b;
          encli_1 += 1;
        }
      }

      li += 1;
    }
  }

  // Update k-mer count.
  kmer_count = encli_0 + encli_1;
  ifs.close();

  // Recording end time.
  auto end = chrono::steady_clock::now();
  cout << "-- Done hashing. Now writing. Time so far: "
       << chrono::duration_cast<chrono::seconds>(end - start).count()
       << " seconds." << endl;

  // Allocate new tag array.
  int8_t *new_tag_arr;
  uint64_t new_tag_arr_size = sigs_row_count * partitions * l;

  try {
    new_tag_arr = new int8_t[new_tag_arr_size];
    cout << "-- Done new tag allocation" << endl;
  } catch (std::bad_alloc &ba) {
    std::cerr << "Exception bad_alloc is caught: " << ba.what() << endl;
  }

  // Unitialize new tag array to -1.
  for (uint64_t j = 0; j < new_tag_arr_size; j++) {
    new_tag_arr[j] = -1;
  }

  // Allocate new enc_ind_array.
  uint8_t *new_encid_arr;
  uint64_t increment = 8;
  uint64_t new_encid_arr_size = ceil(sigs_arr_size / increment);

  try {
    new_encid_arr = new uint8_t[new_encid_arr_size];
    cout << "-- Done new enc id allocation" << endl;
  } catch (std::bad_alloc &ba) {
    std::cerr << "Exception bad_alloc is caught: " << ba.what() << endl;
  }

  // Sort sigs and tag arrays.
  // Sort uint32_t array b[] according to the order defined by a[].
  for (int i = 0; i < l; i++) {
    for (uint64_t j = 0; j < sigs_row_count; j++) {
      if (sigs_indicator_arr[sigs_row_count * i + j] > 0) {
        vector<int8_t> tag_col_count(partitions, 0);
        uint8_t n = sigs_indicator_arr[sigs_row_count * i + j];

        vector<tuple<int8_t, uint32_t, uint8_t>> pairt;

        // Storing the respective array elements in pairs.
        for (uint64_t k = 0; k < n; k++) {
          uint64_t tmp_idx =
              (sigs_row_count * sigs_col_count * partitions * i) +
              (j * sigs_col_count * partitions) + k;
          int8_t tag_val = tag_arr[tmp_idx];
          uint32_t sig_val = sigs_arr[tmp_idx];
          uint8_t encoding_id = enc_arr_id[tmp_idx];

          pairt.push_back(make_tuple(tag_val, sig_val, encoding_id));
          tag_col_count[tag_val] += 1;
        }

        sort(pairt.begin(), pairt.end());

        // Update new_tag_array.
        for (uint64_t r = 0; r < partitions; r++) {
          uint64_t tmp_idx =
              (sigs_row_count * partitions * i) + (j * partitions) + r;
          if (r == 0) {
            new_tag_arr[tmp_idx] = tag_col_count[r];
          } else {
            new_tag_arr[tmp_idx] = new_tag_arr[tmp_idx - 1] + tag_col_count[r];
          }
        }

        // Modifying original sig array.
        for (uint64_t k = 0; k < n; k++) {
          uint64_t tmp_idx =
              (sigs_row_count * sigs_col_count * partitions * i) +
              (j * sigs_col_count * partitions) + k;
          sigs_arr[tmp_idx] = get<1>(pairt[k]);
          enc_arr_id[tmp_idx] = get<2>(pairt[k]);
        }
      }
    }
  }

  // Populate new encid array.
  uint64_t j = 0;
  for (uint64_t i = 0; i < enc_arr_id_size; i += increment) {
    uint8_t combo_id = 0;
    for (uint64_t x = 0; x < increment; x++) {
      combo_id = combo_id << 1;
      combo_id = combo_id | enc_arr_id[i + x];
    }
    new_encid_arr[j] = combo_id;
    j += 1;
  }

  string map_meta = "meta";
  string path = output_library_dir + "/" + map_meta;

  FILE *wfmeta;
  wfmeta = fopen(path.c_str(), "wb");
  if (!wfmeta) {
    cout << "Cannot open file!" << endl;
    return 1;
  }

  // Write metadata file.
  // Write p, l, h values.
  fwrite(&p, sizeof(uint64_t), 1, wfmeta);
  fwrite(&l, sizeof(uint64_t), 1, wfmeta);
  fwrite(&h, sizeof(uint64_t), 1, wfmeta);

  fwrite(&sigs_arr_size, sizeof(uint64_t), 1, wfmeta);
  fwrite(&new_tag_arr_size, sizeof(uint64_t), 1, wfmeta);

  // Write k-mer count.
  fwrite(&kmer_count, sizeof(uint64_t), 1, wfmeta);
  cout << "k-mer count " << kmer_count << endl;
  fwrite(&encli_0, sizeof(uint64_t), 1, wfmeta);
  cout << "k-mer count array 0 " << encli_0 << endl;
  fwrite(&encli_1, sizeof(uint64_t), 1, wfmeta);
  cout << "k-mer count array 1 " << encli_1 << endl;

  fwrite(&new_encid_arr_size, sizeof(uint64_t), 1, wfmeta);

  // Write mask array and output real count k-mers included in DB.
  for (int i = 0; i < l; i++) {
    cout << "k-mers included l = " << i << " >> " << included_kmers_counter[i]
         << endl;

    int size = shifts[i].size();
    fwrite(&size, sizeof(int8_t), 1, wfmeta);

    for (int s = 0; s < size; s++) {
      fwrite(&shifts[i][s], sizeof(int8_t), 1, wfmeta);
    }
    for (int s = 0; s < size; s++) {
      fwrite(&grab_bits[i][s], sizeof(int8_t), 1, wfmeta);
    }
  }

  // Write sig columns counts and  sig row count.
  uint64_t sig_col_cnt = sigs_col_count;
  fwrite(&sig_col_cnt, sizeof(uint64_t), 1, wfmeta);
  fwrite(&sigs_row_count, sizeof(uint64_t), 1, wfmeta);

  // Write tag size, partition count, tag mask and big_sig_mask.
  fwrite(&t, sizeof(uint64_t), 1, wfmeta);
  fwrite(&partitions, sizeof(uint64_t), 1, wfmeta);
  fwrite(&tag_mask, sizeof(uint64_t), 1, wfmeta);
  fwrite(&big_sig_mask, sizeof(uint64_t), 1, wfmeta);

  // Write file chunk information.
  uint8_t sigf_chunks = SIGF_CHUNKS;
  uint8_t tagf_chunks = TAGF_CHUNKS;
  uint8_t encf_chunks = ENCF_CHUNKS;
  fwrite(&sigf_chunks, sizeof(uint8_t), 1, wfmeta);
  fwrite(&tagf_chunks, sizeof(uint8_t), 1, wfmeta);
  fwrite(&encf_chunks, sizeof(uint8_t), 1, wfmeta);

  // Counter variables for sigs, indicator and tags arrays to output counts.
  uint64_t total_sigs_written = 0;
  uint64_t total_tags_written = 0;
  uint64_t total_members_written_0 = 0;
  uint64_t total_members_written_1 = 0;
  uint64_t total_encid_written = 0;

  // Compute members in file chunks.
  vector<uint64_t> sig_chunk_counts(SIGF_CHUNKS, 0);
  vector<uint64_t> tag_chunk_counts(TAGF_CHUNKS, 0);
  vector<uint64_t> enc_chunk_counts_0(ENCF_CHUNKS, 0);
  vector<uint64_t> enc_chunk_counts_1(ENCF_CHUNKS, 0);
  vector<uint64_t> encid_chunk_counts(SIGF_CHUNKS, 0);

  for (int m = 0; m < SIGF_CHUNKS - 1; m++) {
    sig_chunk_counts[m] = round(sigs_arr_size / SIGF_CHUNKS);
    total_sigs_written += sig_chunk_counts[m];
  }
  sig_chunk_counts[SIGF_CHUNKS - 1] = sigs_arr_size - total_sigs_written;

  for (int m = 0; m < TAGF_CHUNKS - 1; m++) {
    tag_chunk_counts[m] = round(new_tag_arr_size / TAGF_CHUNKS);
    total_tags_written += tag_chunk_counts[m];
  }
  tag_chunk_counts[TAGF_CHUNKS - 1] = new_tag_arr_size - total_tags_written;

  for (int m = 0; m < ENCF_CHUNKS - 1; m++) {
    enc_chunk_counts_0[m] = round(encli_0 / ENCF_CHUNKS);
    total_members_written_0 += enc_chunk_counts_0[m];
  }
  enc_chunk_counts_0[ENCF_CHUNKS - 1] = encli_0 - total_members_written_0;

  for (int m = 0; m < ENCF_CHUNKS - 1; m++) {
    enc_chunk_counts_1[m] = round(encli_1 / ENCF_CHUNKS);
    total_members_written_1 += enc_chunk_counts_1[m];
  }
  enc_chunk_counts_1[ENCF_CHUNKS - 1] = encli_1 - total_members_written_1;

  for (int m = 0; m < SIGF_CHUNKS - 1; m++) {
    encid_chunk_counts[m] = round(new_encid_arr_size / SIGF_CHUNKS);
    total_encid_written += encid_chunk_counts[m];
  }
  encid_chunk_counts[SIGF_CHUNKS - 1] =
      new_encid_arr_size - total_encid_written;

  total_sigs_written = 0;
  total_tags_written = 0;
  total_members_written_0 = 0;
  total_members_written_1 = 0;
  total_encid_written = 0;

  // Check total counts and write to file.
  for (int m = 0; m < SIGF_CHUNKS; m++) {
    fwrite(&sig_chunk_counts[m], sizeof(uint64_t), 1, wfmeta);

    // Open sig file.
    string map_sig = "sig" + to_string(m);
    path = output_library_dir + "/" + map_sig;

    FILE *wf;
    wf = fopen(path.c_str(), "wb");

    if (!wf) {
      cout << "Cannot open file!" << endl;
      return 1;
    }

    // Write sigs.
    total_sigs_written += fwrite(sigs_arr + total_sigs_written,
                                 sizeof(uint32_t), sig_chunk_counts[m], wf);
    fclose(wf);
  }

  // Check total counts and write to file.
  for (int m = 0; m < TAGF_CHUNKS; m++) {
    fwrite(&tag_chunk_counts[m], sizeof(uint64_t), 1, wfmeta);

    // Open tag file.
    string map_tag = "tag" + to_string(m);
    path = output_library_dir + "/" + map_tag;

    FILE *wftag;
    wftag = fopen(path.c_str(), "wb");
    if (!wftag) {
      cout << "Cannot open file!" << endl;
      return 1;
    }

    // Write tags.
    total_tags_written += fwrite(new_tag_arr + total_tags_written,
                                 sizeof(int8_t), tag_chunk_counts[m], wftag);
    fclose(wftag);
  }

  // Check total counts and write to file.
  for (int m = 0; m < ENCF_CHUNKS; m++) {
    fwrite(&enc_chunk_counts_0[m], sizeof(uint64_t), 1, wfmeta);

    // Open encoding file.
    string map_enc = "enc" + to_string(m);
    path = output_library_dir + "/" + map_enc;

    FILE *wfenc;
    wfenc = fopen(path.c_str(), "wb");
    if (!wfenc) {
      cout << "Cannot open file!" << endl;
      return 1;
    }

    // Write encoding array.
    total_members_written_0 +=
        fwrite(encode_arr_0 + total_members_written_0, sizeof(uint64_t),
               enc_chunk_counts_0[m], wfenc);
    fclose(wfenc);
  }

  // Check total counts and write to file.
  for (int m = 0; m < ENCF_CHUNKS; m++) {
    fwrite(&enc_chunk_counts_1[m], sizeof(uint64_t), 1, wfmeta);

    // Open encoding file.
    string map_ence = "ence" + to_string(m);
    path = output_library_dir + "/" + map_ence;

    FILE *wfenc;
    wfenc = fopen(path.c_str(), "wb");
    if (!wfenc) {
      cout << "Cannot open file!" << endl;
      return 1;
    }

    // Write encoding array.
    total_members_written_1 +=
        fwrite(encode_arr_1 + total_members_written_1, sizeof(uint64_t),
               enc_chunk_counts_1[m], wfenc);
    fclose(wfenc);
  }

  // Check total encid and write to file.
  for (int m = 0; m < SIGF_CHUNKS; m++) {
    fwrite(&encid_chunk_counts[m], sizeof(uint64_t), 1, wfmeta);

    // Open sig file.
    string map_encid = "encid" + to_string(m);
    path = output_library_dir + "/" + map_encid;

    FILE *wf;
    wf = fopen(path.c_str(), "wb");
    if (!wf) {
      cout << "Cannot open file!" << endl;
      return 1;
    }

    // Write sigs.
    total_encid_written += fwrite(new_encid_arr + total_encid_written,
                                  sizeof(uint8_t), encid_chunk_counts[m], wf);
    fclose(wf);
  }

  cout << "Sigs written : " << total_sigs_written << endl;
  cout << "Tags written : " << total_tags_written << endl;
  cout << "Enc id written : " << total_encid_written << endl;
  cout << "Encodings written to array 0 : " << total_members_written_0 << endl;
  cout << "Encodings written to array 1 : " << total_members_written_1 << endl;
  cout << "Encodings written : "
       << total_members_written_0 + total_members_written_1 << endl;

  fclose(wfmeta);

  end = chrono::steady_clock::now();
  cout << "-- Done making map. Time so far: "
       << chrono::duration_cast<chrono::seconds>(end - start).count()
       << " seconds" << endl;

  // Output map information.
  cout << "Columns = " << sigs_col_count << endl;
  cout << "Partitions = " << unsigned(partitions) << endl;
  cout << "Sigs row count = " << sigs_row_count << endl;
  cout << "Tag size = " << unsigned(t) << endl;
  cout << "Tag mask = " << unsigned(tag_mask) << endl;
  cout << "Big sig mask = " << big_sig_mask << endl;

  // Compute the usage of rows in signature matrix to determine how many
  // positions in each row are used.

  // Row count vector.
  vector<uint64_t> sig_row_count_vec(sigs_col_count * partitions + 1, 0);

  // Traverse indicator array to compute filled positions.
  for (int i = 0; i < l; i++) {
    for (uint64_t r = 0; r < sigs_row_count; r++) {
      sig_row_count_vec[sigs_indicator_arr[sigs_row_count * i + r]] += 1;
    }
    for (int s = 0; s < sigs_col_count * partitions + 1; s++) {
      cout << "l = " << l << " -- Count of rows with positions filled " << s
           << " : " << sig_row_count_vec[s];
      cout << " (" << std::fixed << std::setprecision(6)
           << (double)sig_row_count_vec[s] / sigs_row_count << ")" << endl;
    }
    cout << "l = " << i << " -- Total populated rows : "
         << std::accumulate(sig_row_count_vec.begin(), sig_row_count_vec.end(),
                            decltype(sig_row_count_vec)::value_type(0))
         << endl;
    // Reset counts to 0.
    fill(sig_row_count_vec.begin(), sig_row_count_vec.end(), 0);
  }

  // Do not forget to delete when it is done.
  delete[] sigs_arr;
  delete[] sigs_indicator_arr;
  delete[] tag_arr;
  delete[] encode_arr_0;
  delete[] encode_arr_1;
  delete[] enc_arr_id;
  delete[] new_tag_arr;

  end = chrono::steady_clock::now();
  cout << "-- Done writing. Time so far: "
       << chrono::duration_cast<chrono::seconds>(end - start).count()
       << " seconds" << endl;

  return 0;
}

// Function definitions.

void encodekmer(const char *s, uint64_t &b, uint64_t &b_sig) {
  for (int i = 0; i < int(KMER_LENGTH); i++) {
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

uint64_t file_read(istream &is, vector<char> &buff) {
  is.read(&buff[0], buff.size());
  return is.gcount();
}

uint64_t count_lines(const vector<char> &buff, int size) {
  uint64_t newlines = 0;
  const char *p = &buff[0];
  for (int i = 0; i < size; i++) {
    if (p[i] == '>') {
      newlines++;
    }
  }
  return newlines;
}
