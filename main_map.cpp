#include <algorithm>
#include <bitset>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <new>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define VERSION 17.1
#define KMER_LENGTH 32

#define SIGS_COL_COUNT 7

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

  char *input_fasta_file = NULL;
  char *output_library_dir = NULL;
  uint64_t p_value = 3;

  int cf_tmp;
  opterr = 0;

  while (1) {
    static struct option long_options[] = {
        {"input-fasta-file", 1, 0, 'i'},
        {"output-library-dir", 1, 0, 'o'},
        {"p-value", 1, 0, 'p'},
        {0, 0, 0, 0},
    };

    int option_index = 0;
    cf_tmp = getopt_long(argc, argv, "i:o:p:", long_options, &option_index);

    if ((optarg != NULL) && (*optarg == '-')) {
      cf_tmp = ':';
    }

    if (cf_tmp == -1)
      break;
    switch (cf_tmp) {
    case 'i':
      input_fasta_file = optarg;
      break;
    case 'p':
      p_value = atoi(optarg); // Default value is 3.
      break;
    case 'o':
      output_library_dir = optarg;
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

  if (argc <= 7) {
    printf(
        "-- Arguments supplied are; \nInput file : %s\nOutput directory : %s\n",
        input_fasta_file, output_library_dir);

    const int SZ = 1024 * 1024;
    vector<char> buff(SZ);

    // Open input file.
    string input_fin = input_fasta_file;
    ifstream ifs(input_fin);

    if (!ifs.is_open()) {
      std::cout << "Cannot open file!" << endl;
      exit(0);
    }

    while (uint64_t cc = file_read(ifs, buff)) {
      kmer_count += count_lines(buff, cc);
    }

    ifs.close();

  } else if (argc > 7) {
    printf("Too many arguments are supplied.\n");
    exit(0);
  }

  string input_fin = input_fasta_file;
  ifstream fin(input_fin);
  if (!fin.is_open()) {
    std::cout << "Cannot open file!" << endl;
    exit(0);
  }

  uint8_t k = KMER_LENGTH;
  uint64_t p = p_value;
  uint64_t L = 2;
  /* float alpha = 0.95; // This is just to determine correct h value. */
  uint64_t h = 15;
  /* h = round(log(1 - (pow((1 - alpha), (1 / float(L))))) / */
  /*           log(1 - float(p) / float(k))); */

  cout << "k-mer count = " << kmer_count << endl;
  cout << "k = " << k << '\n';
  cout << "p = " << p << '\n';
  cout << "L = " << L << '\n';
  /* cout << "alpha = " << alpha << '\n'; */
  cout << "Using h = " << h << '\n';

  uint8_t tag_size = 2; // Tag size in bits.
  uint64_t partitions = (pow(2, tag_size));
  uint64_t sigs_row_count = (pow(2, (2 * h) - tag_size));
  uint64_t sigs_col_count = SIGS_COL_COUNT;

  // Allocate signature array.
  uint32_t *sigs_arr;
  uint64_t sigs_arr_size = sigs_row_count * sigs_col_count * partitions * L;

  try {
    sigs_arr = new uint32_t[sigs_arr_size];
    cout << "-- Done sigs allocation" << endl;
  } catch (std::bad_alloc &ba) {
    std::cerr << "Exception bad_alloc is caught: " << ba.what() << endl;
  }

  // Allocate sigs indicator array.
  uint8_t *sigs_indicator_arr;
  uint64_t sigs_indicator_arr_size = sigs_row_count * L;

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
  uint64_t tag_arr_size = sigs_row_count * sigs_col_count * partitions * L;

  try {
    tag_arr = new int8_t[tag_arr_size];
    cout << "-- Done tag allocation" << endl;
  } catch (std::bad_alloc &ba) {
    std::cerr << "Exception bad_alloc is caught: " << ba.what() << endl;
  }

  /* Initialize tag array to -1. */
  // for (uint64_t j = 0; j < sigs_row_count * sigs_col_count * partitions * L;
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

  // Write map to file.
  string map_name = output_library_dir;

  // Create a directory (Linux dependent functionality).
  const int dir_err =
      mkdir(map_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  if (-1 == dir_err) {
    cerr << "Cannot create directory! Error: " << strerror(errno) << endl;
    exit(1);
  }

  // Generate mask.
  uint64_t l, m, n;
  uint64_t included_kmers_counter[L];
  vector<vector<int>> positions(L, vector<int>(h));

  srand(time(NULL));

  for (l = 0; l < L; l++) {
    vector<int> rand_num;

    // Initialize lost k-mers array.
    included_kmers_counter[l] = 0;

    for (m = 0; m < h; m++) {
      n = rand() % k;
      if (count(rand_num.begin(), rand_num.end(), n)) {
        m -= 1;
      } else {
        rand_num.push_back(n);
      }
    }

    cout << "Positions for l = " << l << " >> ";
    sort(rand_num.begin(), rand_num.end(), std::greater<int>());
    for (int j = 0; j < (h - 1); j++) {
      positions[l][j] = rand_num[j];
      cout << positions[l][j] << ", ";
    }
    cout << positions[l][h - 1] << endl;

    rand_num.clear();
  }

  vector<vector<int8_t>> shifts;
  vector<vector<int8_t>> grab_bits;

  for (l = 0; l < L; l++) {
    // Construct a vector of integers.
    vector<int8_t> v;
    vector<int8_t> g;
    int lp = 31;
    int jp = 0;

    for (int p = 0; p < h; p++) {
      if (p == 0) {
        v.push_back((lp - positions[l][p]) * 2);
        lp = positions[l][p];
        jp += 2;

      } else if (positions[l][p - 1] - positions[l][p] != 1) {
        v.push_back((lp - positions[l][p]) * 2);
        lp = positions[l][p];
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
  for (int i = 0; i < tag_size; i++) {
    tag_mask += (pow(2, i)); // 111 if tag = 3
  }
  uint64_t big_sig_mask = 0;
  for (int i = 0; i < ((2 * h) - tag_size); i++) {
    big_sig_mask += (pow(2, i));
  }

  while (std::getline(fin, line)) {
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

      for (l = 0; l < L; l++) {
        sig_hash = encodekmer_bits(b_sig, shifts[l], grab_bits[l]);

        // Get first 2 bits of signature (effectively bits 28 - 27) of 32 bit
        // encoding as tag.
        tag = (sig_hash >> ((2 * h) - tag_size)) & tag_mask;

        // Get last 26 bits of signature (effectively bits 26 - 1 ) of 32 bit
        // encoding as sigs row number.
        big_sig_hash = sig_hash & big_sig_mask;

        uint64_t tmp_idx =
            (sigs_row_count * sigs_col_count * partitions * l) +
            (big_sig_hash * sigs_col_count * partitions) +
            sigs_indicator_arr[sigs_row_count * l + big_sig_hash];

        // Check is row space is available for forward k-mer.
        if (sigs_indicator_arr[sigs_row_count * l + big_sig_hash] <
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
          sigs_indicator_arr[sigs_row_count * l + big_sig_hash] += 1;

          // Increment k-mer counter.
          included_kmers_counter[l] += 1;

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
  fin.close();

  // Recording end time.
  auto end = chrono::steady_clock::now();
  cout << "-- Done hashing. Now writing. Time so far: "
       << chrono::duration_cast<chrono::seconds>(end - start).count()
       << " seconds." << endl;

  // Allocate new tag array.
  int8_t *new_tag_arr;
  uint64_t new_tag_arr_size = sigs_row_count * partitions * L;

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
  for (l = 0; l < L; l++) {
    for (uint64_t j = 0; j < sigs_row_count; j++) {
      if (sigs_indicator_arr[sigs_row_count * l + j] > 0) {
        vector<int8_t> tag_col_count(partitions, 0);
        uint8_t n = sigs_indicator_arr[sigs_row_count * l + j];

        vector<tuple<int8_t, uint32_t, uint8_t>> pairt;

        // Storing the respective array elements in pairs.
        for (uint64_t i = 0; i < n; i++) {
          uint64_t tmp_idx =
              (sigs_row_count * sigs_col_count * partitions * l) +
              (j * sigs_col_count * partitions) + i;
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
              (sigs_row_count * partitions * l) + (j * partitions) + r;
          if (r == 0) {
            new_tag_arr[tmp_idx] = tag_col_count[r];
          } else {
            new_tag_arr[tmp_idx] = new_tag_arr[tmp_idx - 1] + tag_col_count[r];
          }
        }

        // Modifying original sig array.
        for (uint64_t i = 0; i < n; i++) {
          uint64_t tmp_idx =
              (sigs_row_count * sigs_col_count * partitions * l) +
              (j * sigs_col_count * partitions) + i;
          sigs_arr[tmp_idx] = get<1>(pairt[i]);
          enc_arr_id[tmp_idx] = get<2>(pairt[i]);
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

  string map_meta = map_name + "_meta";
  string path = map_name + "/" + map_meta;

  FILE *wfmeta;
  wfmeta = fopen(path.c_str(), "wb");
  if (!wfmeta) {
    cout << "Cannot open file!" << endl;
    return 1;
  }

  // Write metadata file.
  // Write p, L, h values.
  fwrite(&p, sizeof(uint64_t), 1, wfmeta);
  fwrite(&L, sizeof(uint64_t), 1, wfmeta);
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
  for (l = 0; l < L; l++) {
    cout << "k-mers included l = " << l << " >> " << included_kmers_counter[l]
         << endl;

    int size = shifts[l].size();
    fwrite(&size, sizeof(int8_t), 1, wfmeta);

    for (int s = 0; s < size; s++) {
      fwrite(&shifts[l][s], sizeof(int8_t), 1, wfmeta);
    }
    for (int s = 0; s < size; s++) {
      fwrite(&grab_bits[l][s], sizeof(int8_t), 1, wfmeta);
    }
  }

  // Write sig columns counts and  sig row count.
  uint64_t sig_col_cnt = sigs_col_count;
  fwrite(&sig_col_cnt, sizeof(uint64_t), 1, wfmeta);
  fwrite(&sigs_row_count, sizeof(uint64_t), 1, wfmeta);

  // Write tag size, partition count, tag mask and big_sig_mask.
  fwrite(&tag_size, sizeof(uint64_t), 1, wfmeta);
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
    string map_sig = map_name + "_sig" + to_string(m);
    path = map_name + "/" + map_sig;

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
    string map_tag = map_name + "_tag" + to_string(m);
    path = map_name + "/" + map_tag;

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
    string map_enc = map_name + "_enc" + to_string(m);
    path = map_name + "/" + map_enc;

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
    string map_ence = map_name + "_ence" + to_string(m);
    path = map_name + "/" + map_ence;

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
    string map_encid = map_name + "_encid" + to_string(m);
    path = map_name + "/" + map_encid;

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
  cout << "Tag size = " << unsigned(tag_size) << endl;
  cout << "Tag mask = " << unsigned(tag_mask) << endl;
  cout << "Big sig mask = " << big_sig_mask << endl;

  // Compute the usage of rows in signature matrix to determine how many
  // positions in each row are used.

  // Row count vector.
  vector<uint64_t> sig_row_count_vec(sigs_col_count * partitions + 1, 0);

  // Traverse indicator array to compute filled positions.
  for (int l = 0; l < L; l++) {
    for (uint64_t r = 0; r < sigs_row_count; r++) {
      sig_row_count_vec[sigs_indicator_arr[sigs_row_count * l + r]] += 1;
    }
    for (int s = 0; s < sigs_col_count * partitions + 1; s++) {
      cout << "l = " << l << " -- Count of rows with positions filled " << s
           << " : " << sig_row_count_vec[s];
      cout << " (" << std::fixed << std::setprecision(6)
           << (double)sig_row_count_vec[s] / sigs_row_count << ")" << endl;
    }
    cout << "l = " << l << " -- Total populated rows : "
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
