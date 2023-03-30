#include <algorithm>
#include <chrono>
#include <dirent.h>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <math.h>
#include <numeric>
#include <sstream>
#include <string.h> // for strcmp, strlen
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;

#define THREAD_COUNT_OPT 'T'
#define TAXONOMY_LOOKUP_PATH_OPT 'A'
#define KMER_LENGTH 32

namespace TaxonomicInfo {
uint16_t num_ranks = 7;
enum rank { KINGDOM = 1, PHYLUM = 2, CLASS = 3, ORDER = 4, FAMILY = 5, GENUS = 6, SPECIES = 7 };

enum Kingdoms { BACTERIA = 2, ARCHAEA = 2157 };

static const Kingdoms AllKingdoms[] = {BACTERIA, ARCHAEA};

unordered_map<int, string> mapRankName = {{1, "kingdom"}, {2, "phylum"}, {3, "class"},  {4, "order"},
                                          {5, "family"},  {6, "genus"},  {7, "species"}};

} // namespace TaxonomicInfo

struct kmer_match {
  uint16_t dist;
  float vote;
  uint64_t taxID;
};

struct read_info {
  bool isRC;
  string readID;
  vector<kmer_match> match_vector;
};

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

void read_taxonomy_lookup(string filepath, unordered_map<uint64_t, vector<uint64_t>> &taxonomy_lookup) {
  ifstream ftable;
  ftable.open(filepath);
  if (!ftable) {
    cout << "Cannot open file for taxonomic lookup table." << endl;
    exit(1);
  }

  for (string line; getline(ftable, line);) {
    istringstream iss(line);
    string taxIDstr;
    getline(iss, taxIDstr, ' ');
    uint64_t taxID = stoi(taxIDstr);

    vector<uint64_t> lineage;
    string next_taxIDstr;

    while (getline(iss, next_taxIDstr, ',')) {
      lineage.push_back(stoi(next_taxIDstr));
    }
    taxonomy_lookup.insert({taxID, lineage});
  }
}

void read_matches(string filepath, vector<read_info> &all_read_info) {
  ifstream infile(filepath);
  string line;
  string readID;
  uint64_t line_counter = 0;

  uint16_t k = KMER_LENGTH;

  while (getline(infile, line)) {
    istringstream iss(line);

    if ((line_counter % 3) == 0)
      iss >> readID;
    else {
      bool isRC;
      string tmp;
      if ((line_counter % 3) == 1) {
        isRC = false;
        iss >> tmp;
      } else {
        isRC = true;
        iss >> tmp;
      }

      string match_str;
      vector<kmer_match> match_vector;

      while (iss >> match_str) {
        stringstream ss(match_str);
        string taxID_str;
        string dist_str;
        getline(ss, taxID_str, ':');
        getline(ss, dist_str, ':');

        kmer_match curr_match;
        curr_match.dist = stoi(dist_str);
        curr_match.taxID = stoi(taxID_str);
        curr_match.vote = pow((1.0 - curr_match.dist / (float)k), k);
        match_vector.push_back(curr_match);
      }

      read_info curr_read;
      curr_read.readID = readID;
      curr_read.isRC = isRC;
      curr_read.match_vector = match_vector;
      all_read_info.push_back(curr_read);
    }
    line_counter++;
  }
}

unordered_map<uint64_t, float> get_rank_votes(kmer_match amatch,
                                              unordered_map<uint64_t, vector<uint64_t>> &taxonomy_lookup) {
  uint16_t k = KMER_LENGTH;
  unordered_map<uint64_t, float> match_votes;
  for (uint64_t a_taxID : taxonomy_lookup[amatch.taxID]) {
    if ((a_taxID > 0) && (amatch.taxID > 0))
      match_votes[a_taxID] = amatch.vote;
  }
  return match_votes;
}

void aggregate_votes(unordered_map<uint64_t, vector<uint64_t>> taxonomy_lookup,
                     vector<read_info> &all_read_info,
                     unordered_map<uint16_t, unordered_map<uint64_t, float>> &profile_by_rank,
                     int thread_count) {
  int num_reads = all_read_info.size() / 2;
#pragma omp parallel for schedule(dynamic, 1) num_threads(thread_count) shared(taxonomy_lookup, all_read_info)
  for (int rix = 0; rix < num_reads; ++rix) {
    float prev_max_vote = 0.0;
    unordered_map<uint64_t, float> prev_final_votes;
    for (int cx = 0; cx < 2; ++cx) {
      int ix = rix * 2 + cx;
      read_info &curr_read = all_read_info[ix];
      uint64_t rootID = 1;
      pair<uint64_t, float> identity(rootID, 0.0);
      float curr_max_vote = 0.0;
      if (curr_read.match_vector.size() > 0) {
        unordered_map<uint64_t, vector<float>> vote_collector;
        unordered_map<uint64_t, float> curr_final_votes;
        unordered_set<uint64_t> all_taxIDs;
        unordered_map<uint16_t, unordered_set<uint64_t>> taxIDs_by_rank;

        for (auto &curr_match : curr_read.match_vector) {
          unordered_map<uint64_t, float> match_votes;
#pragma omp critical
          { match_votes = get_rank_votes(curr_match, taxonomy_lookup); }
          for (auto &vote : match_votes) {
            vote_collector[vote.first].push_back(vote.second);
            all_taxIDs.insert(vote.first);
#pragma omp critical
            { taxIDs_by_rank[(int)taxonomy_lookup[vote.first].size()].insert(vote.first); }
          }
        }

        vector<uint64_t> taxIDs_vec(all_taxIDs.begin(), all_taxIDs.end());

        for (uint64_t taxID : taxIDs_vec) {
          curr_final_votes[taxID] =
              accumulate(vote_collector[taxID].begin(), vote_collector[taxID].end(), 0.0);
        }
        for (const auto &taxon : TaxonomicInfo::AllKingdoms) {
          curr_max_vote += curr_final_votes[static_cast<int>(taxon)];
        }

        if (ix % 2 == 1) {
          for (uint16_t lvl = TaxonomicInfo::num_ranks; lvl >= 1; --lvl) {
            vector<uint64_t> taxIDs_vec_lvl(taxIDs_by_rank[lvl].begin(), taxIDs_by_rank[lvl].end());
            for (auto &taxID : taxIDs_vec_lvl) {
              if (curr_max_vote > prev_max_vote) {
#pragma omp atomic
                profile_by_rank[lvl][taxID] += curr_final_votes[taxID];
              } else {
#pragma omp atomic
                profile_by_rank[lvl][taxID] += prev_final_votes[taxID];
              }
            }
          }
        } else {
          prev_max_vote = curr_max_vote;
          prev_final_votes = curr_final_votes;
        }
      }
    }
  }
}

void compute_profile(unordered_map<uint16_t, unordered_map<uint64_t, float>> &profile_by_rank,
                     int thread_count) {
#pragma omp parallel for num_threads(min(thread_count, (int)TaxonomicInfo::num_ranks))
  for (uint16_t lvl = TaxonomicInfo::num_ranks; lvl >= 1; --lvl) {
    float total_rank_nvote = 0.0;
    for (auto &kv : profile_by_rank[lvl]) {
      total_rank_nvote += sqrt(kv.second);
    }
    for (auto &kv : profile_by_rank[lvl]) {
      profile_by_rank[lvl][kv.first] = sqrt(kv.second) / total_rank_nvote;
    }
  }
}

void write_profile_to_file(string filepath, string taxonomy_lvl, unordered_map<uint64_t, float> &profile) {
  filepath = filepath + "-" + taxonomy_lvl;
  ofstream outfile(filepath);
  outfile << "TAXONOMY_ID"
          << "\t"
          << "TAXONOMY_LEVEL"
          << "\t"
          << "FRACTION_TOTAL" << endl;
  for (auto &kv : profile) {
    outfile << kv.first << "\t" << taxonomy_lvl << "\t" << to_string(kv.second) << endl;
  }
}

int main(int argc, char *argv[]) {
  uint64_t thread_count = 1;
  char *input_path = NULL;
  string taxonomy_lookup_path;
  string output_predictions_dir = ".";

  int cf_tmp;
  opterr = 0;

  while (1) {
    static struct option long_options[] = {
        {"input-matches-dir", 1, 0, 'i'},
        {"output-predictions-dir", 1, 0, 'o'},
        {"taxonomy-lookup-path", 1, 0, TAXONOMY_LOOKUP_PATH_OPT},
        {"thread-count", 1, 0, THREAD_COUNT_OPT},
        {0, 0, 0, 0},
    };

    int option_index = 0;
    cf_tmp = getopt_long(argc, argv, "i:o:", long_options, &option_index);

    if ((optarg != NULL) && (*optarg == '-')) {
      cf_tmp = ':';
    }

    if (cf_tmp == -1)
      break;
    else if (cf_tmp == TAXONOMY_LOOKUP_PATH_OPT)
      taxonomy_lookup_path = optarg;
    else if (cf_tmp == THREAD_COUNT_OPT)
      thread_count = atoi(optarg); // Default is 1.
    else {
      switch (cf_tmp) {
      case 'i':
        input_path = optarg;
        break;
      case 'o':
        output_predictions_dir = optarg;
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
  cout << "Number of threads is " << thread_count << "." << endl;

  struct stat s_input_path;
  vector<string> match_info_path_list;
  if (stat(input_path, &s_input_path) == 0) {
    if (s_input_path.st_mode & S_IFDIR) {
      match_info_path_list = list_dir(input_path);
    } else if (s_input_path.st_mode & S_IFREG) {
      match_info_path_list.push_back(input_path);
    } else {
      cout << "Filetype in the given query path is not recognized." << endl;
      exit(1);
    }
  } else {
    cout << "Given query path is not valid!" << endl;
    exit(1);
  }

  int total_numberOfReads = 0;
  int total_readTime = 0;
  int total_profilingTime = 0;

  chrono::steady_clock::time_point t1;
  chrono::steady_clock::time_point t2;

  for (string input_path : match_info_path_list) {
    string query_name = input_path.substr(input_path.find_last_of("/") + 1);
    query_name = query_name.substr(0, query_name.find_last_of("."));
    query_name = query_name.substr(query_name.find_first_of("_") + 1);
    string output_path = output_predictions_dir + "/" + "profile_" + query_name;
    cout << "Now processing: " << query_name << endl;

    unordered_map<uint64_t, vector<uint64_t>> taxonomy_lookup;
    read_taxonomy_lookup(taxonomy_lookup_path, taxonomy_lookup);

    vector<read_info> all_read_info;
    t1 = chrono::steady_clock::now();
    read_matches(input_path, all_read_info);
    t2 = chrono::steady_clock::now();

    total_readTime += chrono::duration_cast<chrono::milliseconds>(t2 - t1).count();
    total_numberOfReads += all_read_info.size();

    unordered_map<uint16_t, unordered_map<uint64_t, float>> profile_by_rank;
    t1 = chrono::steady_clock::now();
    aggregate_votes(taxonomy_lookup, all_read_info, profile_by_rank, thread_count);
    compute_profile(profile_by_rank, thread_count);
    t2 = chrono::steady_clock::now();

    total_profilingTime += chrono::duration_cast<chrono::milliseconds>(t2 - t1).count();

    for (uint16_t lvl = TaxonomicInfo::num_ranks; lvl >= 1; --lvl) {
      write_profile_to_file(output_path, TaxonomicInfo::mapRankName[lvl], profile_by_rank[lvl]);
    }
  }
  cout << "Time past (read_matches) = " << total_readTime << "[ms]" << endl;
  cout << "Time past (aggregate_votes + compute_profile) = " << total_profilingTime << "[ms]" << endl;
  cout << "Total number of reads processed: " << total_numberOfReads / 2 << endl;
}
