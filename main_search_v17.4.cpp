// main.cpp -- defining, prototyping, and calling a function

//to compile:
// export PATH=/usr/local/bin:$PATH

// g++-9  main_search_v17.4.cpp -std=c++11 -fopenmp -O3 -o main_search
//on server: g++  main_search_v17.2.cpp -std=c++11 -fopenmp -O3 -o main_search91

// tag is taken from the most significant bits of the kmer
// based on version 14.1 plus adding two encoding arrays
// outputs UCSEQ

//to run:
// ./main_search -i G000307305_nbr_map -c 0 -t 4 -q /Users/admin/CLionProjects/hamming_search_1.0/excluded_fna_fq_downSmpl10M
// ./main_search -i G000399765_nbr_map -c 0 -t 4  -q /Users/admin/CLionProjects/hamming_search_1.0/excluded_fna_fq_downSmpl10M
// ./main_search -i G000016385_nbr_map -c 0 -t 4 -q /Users/admin/CLionProjects/hamming_search_1.0/excluded_fna_fq_downSmpl10M



// brew install libomp

// works with map from main_map_v9.1

//can NOT BE TESTED LOCALLY because of memory allocation


#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iterator>
#include <vector>
#include <map>
#include <string>
#include <cstring>
#include <algorithm>
#include <bitset>
#include <unordered_map>
#include <time.h>
#include <dirent.h>
#include <string.h>
#include <sys/types.h>
#include <chrono>
#include <new>
#include <stdio.h>
#include <omp.h>

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define SL 32
//#define SIGS_COLMN 6


using namespace std;

// prototypes
void encodekmer(const char s[], uint64_t &b, uint64_t &b_sig);
void encodekmer_rev(const char s[], uint64_t &b, uint64_t &b_sig);

void update_kmer(const char *s, uint64_t &b, uint64_t &b_sig);
void update_kmer_rev(const char *s, uint64_t &b, uint64_t &b_sig);

uint64_t encodekmer_bits(uint64_t val, vector<int8_t> shifts, vector<int8_t> bits_to_grab);
uint64_t encodekmer_bits_rev(const char *s, vector<int> pos);

bool hd(uint64_t x, uint64_t y, int m);
vector<string> list_dir(const char *path);
uint8_t get_encid(uint64_t sind, uint8_t enc_arr_id[]);



int main(int argc, char *argv[]) {

    //typedef unsigned long long int  uint64_t;
    auto start = chrono::steady_clock::now();


    char *ivalue = NULL;
    uint64_t cvalue = 0;
    uint64_t thread_count = 4;
    char *qvalue = NULL;
    //int index;
    int cf;

    opterr = 0;


    while ((cf = getopt (argc, argv, "i:c:t:q:")) != -1)
        switch (cf)
        {
            case 'i':
                ivalue = optarg;
                //printf("Option -i has arg: %s\n", optarg);
                break;
            case 'c':
                cvalue = atoi(optarg); // default is 0
                printf("Option -c has arg: %s\n", optarg);
                break;
            case 't':
                thread_count = atoi(optarg); // default is 4
                //printf("Option -t has arg: %s\n", optarg);
                break;
            case 'q':
                qvalue = optarg;
                //printf("Option -fq has arg: %s\n", optarg);
                break;
            case '?':
                if (optopt == 'i')
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);

                else if (optopt == 'q')
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);

                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);

                else
                    fprintf (stderr,
                             "Unknown option character `\\x%x'.\n",
                             optopt);
                return 1;
            default:
                abort ();
        }






    if (argc <= 9)
    {
        //printf("-- Arguments supplied are \n");
        printf ("-- Arguments supplied are \nMap %s \nQuery %s\nThreads %d\n", ivalue, qvalue, thread_count);
        //printf("Filename %s\n", argv[2]);
        //printf("p = %d\n", atoi(argv[4]));
    }
    else if (argc > 9) {
        printf("Too many arguments supplied.\n");
        exit(0);
    }

    //read map

    string map_name = ivalue;

    string map_meta = map_name + "_meta";
    string path = map_name + "/" + map_meta;


    FILE *fmeta = fopen(path.c_str(), "rb");
    if (! fmeta) {
        cout << "Cannot open file!" << endl;
        return 1;
    }



    // set confidence threshold, should be at 0
    uint64_t c = cvalue;
    cout << "cvalue " << uint64_t(cvalue) << "\n";


    // read parameters from input file
    uint64_t p;
    fread(&p, sizeof(uint64_t), 1, fmeta);

    uint64_t L;
    fread(&L, sizeof(uint64_t), 1, fmeta);

    float alpha;
    fread(&alpha, sizeof(float), 1, fmeta);

    uint64_t K;
    fread(&K, sizeof(uint64_t), 1, fmeta);

    uint64_t sigs_arr_size;
    fread(&sigs_arr_size, sizeof(uint64_t), 1, fmeta);

    uint64_t new_tag_arr_size;
    fread(&new_tag_arr_size, sizeof(uint64_t), 1, fmeta);

    uint64_t kmer_count;
    fread(&kmer_count, sizeof(uint64_t), 1, fmeta);

    uint64_t encli_0;
    fread(&encli_0, sizeof(uint64_t), 1, fmeta);
    uint64_t encli_1;
    fread(&encli_1, sizeof(uint64_t), 1, fmeta);

    uint64_t enc_arr_id_size;
    fread(&enc_arr_id_size, sizeof(uint64_t), 1, fmeta);


    string line;
    uint64_t l;


    //cout << "Map " << map_name << endl;
    cout << "k-mer count = " << kmer_count << endl;
    cout << "kmer count array 0  = " << encli_0 << endl;
    cout << "kmer count array 1  = " << encli_1 << endl;
    cout << "Sig array size = " << sigs_arr_size << endl;
    cout << "Tag array size  = " << new_tag_arr_size << endl;
    cout << "Enc id array size  = " << enc_arr_id_size << endl;
    cout << "p = " << p << '\n';
    cout << "SL = " << SL << '\n';
    cout << "L = " << L << '\n';
    cout << "alpha = " << alpha << '\n';
    cout << "Using K = " << K << '\n';


    // read mask
    // vector<vector<int> > positions(L, vector<int>(K));

    vector<vector<int8_t>> shifts;
    vector<vector<int8_t>> grab_bits;
    int vec_size = 0;


    for (l = 0; l < L; l++)
    {

        vector<int8_t> v;
        vector<int8_t> g;
        int8_t val_read;

        fread(&vec_size, sizeof(int8_t), 1, fmeta);
        //cout << vec_size << endl;

        for (int s = 0; s < vec_size; s++)
        {
            fread(&val_read, sizeof(int8_t), 1, fmeta);
            //cout << "read " << signed(val_read) << endl;
            v.push_back(val_read);
        }
        shifts.push_back(v);
        v.clear();



        for (int s = 0; s < vec_size; s++)
        {
            fread(&val_read, sizeof(int8_t), 1, fmeta);
            g.push_back(val_read);
        }
        grab_bits.push_back(g);
        g.clear();

    }


//    int q  = 0;
//    while (shifts[0][q] !=-1)
//    {
//        cout << "shifts " << signed(shifts[0][q]) << endl;
//        cout << "grab_bits " << signed(grab_bits[0][q]) << endl;
//        q++;
//    }



    // read sig columns counts and  sig row count
    uint64_t SIGS_COLMN;
    fread(&SIGS_COLMN, sizeof(uint64_t), 1, fmeta);

    uint64_t sigs_row_count;
    fread(&sigs_row_count, sizeof(uint64_t), 1, fmeta);
    //cout << sigs_row_count << endl;


    // read tag size, partition count, tag mask and big_sig_mask
    uint64_t tag_size;
    fread(&tag_size, sizeof(uint64_t), 1, fmeta);
    uint64_t partitions;
    fread(&partitions, sizeof(uint64_t), 1, fmeta);
    uint64_t tag_mask;
    fread(&tag_mask, sizeof(uint64_t), 1, fmeta);
    uint64_t big_sig_mask;
    fread(&big_sig_mask, sizeof(uint64_t), 1, fmeta);

    //cout << unsigned(SIGS_COLMN) << " " << sigs_row_count << " " << unsigned(tag_size) << " " << unsigned(partitions) << " " << unsigned(tag_mask) << " " << big_sig_mask << endl;


    // output map information
    cout << "Columns = " << unsigned(SIGS_COLMN) << endl;
    cout << "Partitions = " << unsigned(partitions) << endl;
    cout << "Sigs row count = " << sigs_row_count << '\n';
    cout << "Tag size = " << unsigned(tag_size) << '\n';
    cout << "Tag mask = " << unsigned(tag_mask) << '\n';
    cout << "Big sig mask = " << big_sig_mask << '\n';


    // read file chunk information
    uint8_t sigf_chunks;
    fread(&sigf_chunks, sizeof(uint8_t), 1, fmeta);
    uint8_t tagf_chunks;
    fread(&tagf_chunks, sizeof(uint8_t), 1, fmeta);
    uint8_t encf_chunks;
    fread(&encf_chunks, sizeof(uint8_t), 1, fmeta);



    // read members in file chunks
    vector<uint64_t> sig_chunk_cnts;
    vector<uint64_t> tag_chunk_cnts;
    vector<uint64_t> enc_chunk_cnts;
    vector<uint64_t> enc_chunk_cnts_1;
    vector<uint64_t> encid_chunk_cnts;

    vector<uint64_t> sig_chunk_strt;
    vector<uint64_t> tag_chunk_strt;
    vector<uint64_t> enc_chunk_strt;
    vector<uint64_t> enc_chunk_strt_1;
    vector<uint64_t> encid_chunk_strt;



    uint64_t read_counts;
    uint64_t sum_read_counts = 0;

    for (int m = 0; m < sigf_chunks; m++){
        fread(&read_counts, sizeof(uint64_t), 1, fmeta);
        sig_chunk_cnts.push_back(read_counts);
        sig_chunk_strt.push_back(sum_read_counts);
        sum_read_counts += read_counts;
    }


    sum_read_counts = 0;

    for (int m = 0; m < tagf_chunks; m++){
        fread(&read_counts, sizeof(uint64_t), 1, fmeta);
        tag_chunk_cnts.push_back(read_counts);
        tag_chunk_strt.push_back(sum_read_counts);
        sum_read_counts += read_counts;
    }


    sum_read_counts = 0;

    for (int m = 0; m < encf_chunks; m++){
        fread(&read_counts, sizeof(uint64_t), 1, fmeta);
        enc_chunk_cnts.push_back(read_counts);
        enc_chunk_strt.push_back(sum_read_counts);
        sum_read_counts += read_counts;
    }


    sum_read_counts = 0;

    for (int m = 0; m < encf_chunks; m++){
        fread(&read_counts, sizeof(uint64_t), 1, fmeta);
        enc_chunk_cnts_1.push_back(read_counts);
        enc_chunk_strt_1.push_back(sum_read_counts);
        sum_read_counts += read_counts;
    }

    sum_read_counts = 0;

    for (int m = 0; m < sigf_chunks; m++){
        fread(&read_counts, sizeof(uint64_t), 1, fmeta);
        encid_chunk_cnts.push_back(read_counts);
        encid_chunk_strt.push_back(sum_read_counts);
        sum_read_counts += read_counts;
    }

    fclose(fmeta);


    // allocate signature array

    uint32_t * sigs_arr;

    try
    {
        sigs_arr = new uint32_t [sigs_arr_size];
        cout << "-- Done sigs allocation " << endl;
    }
    catch (std::bad_alloc& ba)
    {
        std::cerr << "bad_alloc caught: " << ba.what() << '\n';
    }



    // allocate sigs indicator array

//    uint8_t * sigs_indicator_arr;
//
//    try
//    {
//        sigs_indicator_arr = new uint8_t [sigs_row_count*L];
//        cout << "-- Done indicator allocation " << endl;
//    }
//    catch (std::bad_alloc& ba)
//    {
//        std::cerr << "bad_alloc caught: " << ba.what() << '\n';
//    }



    // allocate tag array

    int8_t * tag_arr;

    try
    {
        tag_arr = new int8_t [new_tag_arr_size];
        cout << "-- Done tag allocation " << endl;
    }
    catch (std::bad_alloc& ba)
    {
        std::cerr << "bad_alloc caught: " << ba.what() << '\n';
    }



    // allocate and read encoding array

    uint64_t *encode_arr;
    uint64_t *encode_arr_1;

    try
    {
        encode_arr = new uint64_t[encli_0];
        cout << "-- Done encoding array 0 allocation " << endl;
    }
    catch (std::bad_alloc &ba)
    {
        std::cerr << "bad_alloc caught: " << ba.what() << '\n';
    }

    try
    {
        encode_arr_1 = new uint64_t[encli_1];
        cout << "-- Done encoding array 1 allocation " << endl;
    }
    catch (std::bad_alloc &ba)
    {
        std::cerr << "bad_alloc caught: " << ba.what() << '\n';
    }


    // allocate enc array id array
    uint8_t * enc_arr_id;

    try
    {
        enc_arr_id = new uint8_t [enc_arr_id_size];
        cout << "-- Done enc id allocation " << endl;
    }
    catch (std::bad_alloc& ba)
    {
        std::cerr << "bad_alloc caught: " << ba.what() << '\n';
    }

    uint64_t total_sigs_read = 0;
    //uint64_t total_indicators_read = 0;
    uint64_t total_tags_read = 0;
    uint64_t num_pairs = 0;
    uint64_t num_pairs_1 = 0;
    uint64_t total_encid_read = 0;



//    FILE *fmeta = fopen(map_meta.c_str(), "rb");
//    FILE *ftag = fopen(map_tag.c_str(), "rb");
//    FILE *fenc = fopen(map_enc.c_str(), "rb");


    // read files in parallel

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

    // read signatures array
    #pragma omp parallel num_threads(thread_count) shared (sigs_arr, total_sigs_read)
    {
        #pragma omp for
        for (int m = 0; m < sigf_chunks; m++) {
            // open sig file
//            string map_sig = map_name + "_sig" + to_string(m);
//            path = map_name + "/" + map_sig;

            FILE *f;
            f = fopen(str_map_sig[m].c_str(), "rb");

            if (!f) {
                cout << "Cannot open file1!" << endl;
                exit(0);
            }

            // read sigs
            temp_count = fread(sigs_arr + sig_chunk_strt[m], sizeof(uint32_t), sig_chunk_cnts[m], f);

            #pragma omp atomic
            total_sigs_read += temp_count;

            fclose(f);
        }
    }


    // read tag array
    #pragma omp parallel num_threads(thread_count) shared (tag_arr, total_tags_read)
    {
        #pragma omp for
        for (int m = 0; m < tagf_chunks; m++) {
            // open tag file
//            string map_tag = map_name + "_tag" + to_string(m);
//            path = map_name + "/" + map_tag;

            FILE *ftag;
            ftag = fopen(str_map_tag[m].c_str(), "rb");

            if (!ftag) {
                cout << "Cannot open file!2" << endl;
                exit(0);
            }

            // read tags
            temp_count =fread(tag_arr + tag_chunk_strt[m], sizeof(int8_t), tag_chunk_cnts[m], ftag);


            #pragma omp atomic
            total_tags_read += temp_count;

            fclose(ftag);
        }
    }


    // read encoding array
    #pragma omp parallel num_threads(thread_count) shared (encode_arr, num_pairs)
    {
        #pragma omp for
        for (int m = 0; m < encf_chunks; m++) {
            // open encode file
//            string map_enc = map_name + "_enc" + to_string(m);
//            path = map_name + "/" + map_enc;

            FILE *fenc;
            fenc = fopen(str_map_enc[m].c_str(), "rb");

            if (!fenc) {
                cout << "Cannot open file3!" << endl;
                exit(0);
            }

            // read encodings
            temp_count = fread(encode_arr + enc_chunk_strt[m], sizeof(uint64_t), enc_chunk_cnts[m], fenc);

            #pragma omp atomic
            num_pairs += temp_count;

            fclose(fenc);
        }
    }


    // read encoding array
    #pragma omp parallel num_threads(thread_count) shared (encode_arr_1, num_pairs_1)
    {
    #pragma omp for
        for (int m = 0; m < encf_chunks; m++) {
            // open encode file

            FILE *fence;
            fence = fopen(str_map_ence[m].c_str(), "rb");

            if (!fence) {
                cout << "Cannot open file3!" << endl;
                exit(0);
            }

            // read encodings
            temp_count = fread(encode_arr_1 + enc_chunk_strt_1[m], sizeof(uint64_t), enc_chunk_cnts_1[m], fence);

            #pragma omp atomic
            num_pairs_1 += temp_count;

            fclose(fence);
        }
    }

    // read encid array
    #pragma omp parallel num_threads(thread_count) shared (enc_arr_id, total_encid_read)
    {
    #pragma omp for
        for (int m = 0; m < sigf_chunks; m++) {
            // open sig file
//            string map_sig = map_name + "_sig" + to_string(m);
//            path = map_name + "/" + map_sig;

            FILE *fencid;
            fencid = fopen(str_map_encid[m].c_str(), "rb");


            if (!fencid) {
                cout << "Cannot open file5!" << endl;
                exit(0);
            }

            // read sigs
            temp_count = fread(enc_arr_id + encid_chunk_strt[m], sizeof(uint8_t), encid_chunk_cnts[m], fencid);

            #pragma omp atomic
            total_encid_read +=temp_count;

            fclose(fencid);
        }
    }


    //total_indicators_read = fread(sigs_indicator_arr, sizeof(uint8_t), sigs_row_count*L, f);





    cout << "Sigs read  " << ": " << total_sigs_read << endl;
    //cout << "Indicators read  " << ": " << total_indicators_read << endl;
    cout << "Tags read  " << ": " << total_tags_read << endl;
    cout << "Encodings 0 read " << " : " << num_pairs << endl;
    cout << "Encodings 1 read " << " : " << num_pairs_1 << endl;
    cout << "Enc id read " << " : " << total_encid_read << endl;
    //cout << encode_arr[0] << endl;



    auto end = chrono::steady_clock::now();
    cout << "-- Done reading. Now matching. Time so far: "
         << chrono::duration_cast<chrono::seconds>(end - start).count() << " seconds" << endl;



    // read input fastq
    const char *dir = qvalue;


    vector<string> file_list;
    if (dir == NULL) {
        return (1);
    } else {
        file_list = list_dir(dir);
    }

    int file_count = file_list.size();



//    #pragma omp parallel
//        {
//            // Code inside this region runs in parallel.
//            printf("Hello!\n");
//        }

    #pragma omp parallel num_threads(thread_count)
    {
        #pragma omp for schedule(dynamic)
        for (uint64_t f = 0; f < file_count; f++) {


            string input_fq = file_list[f];
            //cout << input_fq << endl;

            stringstream stream;
            stream << std::fixed << std::setprecision(2) << alpha;
            std::string alpha_s = stream.str();

            string input_fq_truct = input_fq.substr(input_fq.find_last_of("/") + 1);

            //string input_map = argv[2];
            //input_map = input_map.substr(0, input_map.find("_"));

            //    string output_fname = "test_k" + to_string(int(SL)) + "_l" + to_string(int(L)) + "_p" +
            //            to_string(int(p))+ "_alpha" + alpha_s + "_" + input_map + "_" +input_fq_truct;

//            string output_fname = "test_k" + to_string(int(SL)) + "_l" + to_string(int(L)) + "_p" +
//                                  to_string(int(p)) + "_alpha" + alpha_s + "_" + input_fq_truct;

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
            uint8_t enc_arr_indict;

            uint64_t kmer_sig;
            int8_t tag;
            uint64_t big_sig_hash;
            uint64_t enc_start;
            uint64_t enc_end;

            while (!ifs.eof()) {
                getline(ifs, line_of_file);


                if (ifs) {

                    if (lines_read % 4 == 0) {
                        name = line_of_file;
                    }

                    if (lines_read % 4 == 1) {
                        int matched = 0;
                        //cout << lines_read << endl;
                        string line_of_file_orig = line_of_file;

                        istringstream iss(line_of_file);

                        while (getline(iss, token, 'N'))
                        {
                            //std::cout << token << std::endl;


                            if (token.length() >= int(SL)){

                                b = 0;
                                b_sig = 0;


                            for (uint64_t i = 0; i < token.length(); i++) {
                                //cout << "i " << i << endl;
                                //cout << "letter " << line_of_file[i]<< endl;
                                //cout << kmer_str << endl;
                                //cout << line_of_file.length() << endl;


                                if (i == 0) {
                                    string kmer_str = token.substr(i, int(SL));
                                    const char *ckmer = kmer_str.c_str();
                                    encodekmer(ckmer, b, b_sig);
                                    i = SL - 1;

                                } else {
                                    string kmer_str = token.substr(i, 1);
                                    //cout << kmer_str << endl;
                                    const char *ckmer = kmer_str.c_str();
                                    update_kmer(ckmer, b, b_sig);
                                }

//                        cout << "b " << std::bitset<64>(b) << endl;
//                        cout << "b_sig " << std::bitset<64>(b_sig) << endl;



                                bool kmerfound = false;

                                for (int64_t funci = 0; funci < L; funci++) {
                                    //cout << "funci " << funci << endl;

                                    kmer_sig = encodekmer_bits(b_sig, shifts[funci], grab_bits[funci]);
                                    //cout << kmer_sig << endl;

                                    // get first 2 bits of signature (effectively bits 28 - 27) of 32 bit encoding as partition numbers
                                    //std::cout << "sig_hash = " << std::bitset<32>(kmer_sig)  << std::endl;

                                    tag = (kmer_sig >> ((2 * K) - tag_size)) & tag_mask;
                                    //std::cout << "tag = " << std::bitset<32>(tag)  << std::endl;
                                    //cout << "tag " << unsigned(tag) << endl;


                                    // get last 26 bits of signature (effectively bits 26 -1 ) of 32 bit encoding as sigs row number
                                    big_sig_hash = kmer_sig & big_sig_mask;
                                    //std::cout << "big_sig = " << std::bitset<32>(big_sig_hash)  << std::endl;




                                    if (tag_arr[sigs_row_count * partitions * funci +
                                                big_sig_hash * partitions + tag] != -1) {
                                        if (tag == 0) {
                                            enc_start = 0;
                                        } else {
                                            enc_start = tag_arr[sigs_row_count * partitions * funci +
                                                                big_sig_hash * partitions + tag - 1];
                                        }
                                        enc_end = tag_arr[sigs_row_count * partitions * funci +
                                                          big_sig_hash * partitions + tag];

                                        //cout << enc_start << " : " << enc_end << endl;



                                        for (uint64_t enc = enc_start; enc < enc_end; enc++) {

                                            // set encoding array id


                                            enc_arr_indict = get_encid(
                                                    sigs_row_count * SIGS_COLMN * partitions * funci +
                                                    big_sig_hash * SIGS_COLMN * partitions + enc,
                                                    enc_arr_id);


                                            if (enc_arr_indict == 0) {
                                                test_enc = encode_arr[sigs_arr[
                                                        sigs_row_count * SIGS_COLMN * partitions * funci +
                                                        big_sig_hash * SIGS_COLMN * partitions + enc]];
                                            } else {
                                                test_enc = encode_arr_1[sigs_arr[
                                                        sigs_row_count * SIGS_COLMN * partitions * funci +
                                                        big_sig_hash * SIGS_COLMN * partitions + enc]];
                                            }

                                            if (hd(b, test_enc, p)) {
                                                matched += 1;
                                                kmerfound = true;
                                                break;
                                            }
                                        }

                                    }


                                    if (kmerfound) {
                                        break;
                                    }
                                }


                                if (matched > c) {
                                    break;
                                }

                            }
                        }

                            if (matched > c) {
                                break;
                            }

                    }



                        // try reverse complement
                        if (matched == 0) {

                            int len = strlen(line_of_file.c_str());
                            char swap;

                            for (int i = 0; i < len / 2; i++) {
                                swap = line_of_file[i];
                                line_of_file[i] = line_of_file[len - i - 1];
                                line_of_file[len - i - 1] = swap;
                            }

                            istringstream iss(line_of_file);

                            while (getline(iss, token, 'N')){
                                //std::cout << token << std::endl;


                                if (token.length() >= int(SL)){


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
                                            //cout << kmer_str << endl;
                                            const char *ckmer = kmer_str.c_str();
                                            update_kmer_rev(ckmer, b, b_sig);
                                        }

                                        //cout << b << endl;


                                        bool kmerfound = false;

                                        for (uint64_t funci = 0; funci < L; funci++) {
                                            //cout << "funci " << funci << endl;

                                            kmer_sig = encodekmer_bits(b_sig, shifts[funci], grab_bits[funci]);
                                            //cout << kmer_sig << endl;
                                            //std::cout << "sig_hash = " << std::bitset<32>(kmer_sig)  << std::endl;

                                            // get first 2 bits of signature (effectively bits 28 - 27) of 32 bit encoding as partition numbers
                                            tag = (kmer_sig >> ((2 * K) - tag_size)) & tag_mask;
                                            //std::cout << "tag = " << std::bitset<8>(tag)  << std::endl;
                                            //cout << "tag " << unsigned(tag) << endl;


                                            // get last 26 bits of signature (effectively bits 26 -1 ) of 32 bit encoding as sigs row number
                                            big_sig_hash = kmer_sig & big_sig_mask;
                                            //cout << big_sig_hash << endl;
                                            //std::cout << "big_sig = " << std::bitset<32>(big_sig_hash)  << std::endl;


                                            if (tag_arr[sigs_row_count * partitions * funci +
                                                        big_sig_hash * partitions + tag] != -1) {
                                                if (tag == 0) {
                                                    enc_start = 0;
                                                } else {
                                                    enc_start = tag_arr[sigs_row_count * partitions * funci +
                                                                        big_sig_hash * partitions + tag - 1];
                                                }
                                                enc_end = tag_arr[sigs_row_count * partitions * funci +
                                                                  big_sig_hash * partitions + tag];


                                                for (uint64_t enc = enc_start; enc < enc_end; enc++) {

                                                    // set encoding array id

                                                    enc_arr_indict = get_encid(
                                                            sigs_row_count * SIGS_COLMN * partitions * funci +
                                                            big_sig_hash * SIGS_COLMN * partitions + enc,
                                                            enc_arr_id);

                                                    if (enc_arr_indict == 0) {
                                                        test_enc = encode_arr[sigs_arr[
                                                                sigs_row_count * SIGS_COLMN * partitions * funci +
                                                                big_sig_hash * SIGS_COLMN * partitions + enc]];
                                                    } else {
                                                        test_enc = encode_arr_1[sigs_arr[
                                                                sigs_row_count * SIGS_COLMN * partitions * funci +
                                                                big_sig_hash * SIGS_COLMN * partitions + enc]];
                                                    }

                                                    if (hd(b, test_enc, p)) {
                                                        matched += 1;
                                                        kmerfound = true;
                                                        break;
                                                    }
                                                }
                                            }


                                            if (kmerfound) {
                                                break;
                                            }
                                        }
                                        if (matched > c) {
                                            break;
                                        }
                                    }
                                }
                                if (matched > c) {
                                    break;
                                }

                            }
                        }


                        if (matched <= c) {
                            //outputFile << ">" << name << " " << matched << endl;
                            outputFile << name << " " << matched << endl;
                            //outputFile << name << endl;
                            //outputFile << line_of_file << endl;
                            outputFile << line_of_file_orig << endl;

                            // output separator and quality
                            getline(ifs, line_of_file);
                            ++lines_read;
                            outputFile << line_of_file << endl;

                            getline(ifs, line_of_file);
                            ++lines_read;
                            outputFile << line_of_file << endl;



                        }
                        else if (matched > c)
                        {
                            reads_matched += 1;
                        }

                    }


//                    if (lines_read % 100000 == 0) {
//                        cout << lines_read << " " << reads_matched << endl;
//                    }

                    ++lines_read;
                } else break;
            }

            #pragma omp critical
            {
                cout << input_fq << " " << lines_read << " " << reads_matched << endl;
            }

//            end = chrono::steady_clock::now();
//            cout << "-- Done matching. Time so far: " << chrono::duration_cast<chrono::seconds>(end - start).count()
//                 << " seconds" << endl;

            ifs.close();
            outputFile.close();
        }
    }

    delete[] sigs_arr; // remember to delete array when I am done
    //delete [] sigs_indicator_arr ;
    delete [] tag_arr ;
    delete[] encode_arr; // remember to delete array when I am done
    delete[] encode_arr_1; // remember to delete array when I am done
    delete[] enc_arr_id;


    end = chrono::steady_clock::now();
    cout << "-- Done matching for all. Time so far: " << chrono::duration_cast<chrono::seconds>(end - start).count()
         << " seconds" << endl;

    return 0;

    }



// function definition

void encodekmer(const char *s, uint64_t &b, uint64_t &b_sig){

    using namespace std;

    for (int i = 0; i <int(SL); i++)
    {
        b = b << 1;
        b_sig = b_sig << 2;

        if (s[i] == 'T'){
            b += 4294967297;
            b_sig += 3;
            //std::cout << "d = " << std::bitset<64>(d)  << std::endl;
        }

        else if (s[i] == 'G'){
            b +=4294967296;
            b_sig +=2;
            //std::cout << "d = " << std::bitset<64>(d)  << std::endl;
        }

        else if (s[i] == 'C'){
            b +=1;
            b_sig +=1;
            //std::cout << "d = " << std::bitset<64>(d)  << std::endl;
        }

        else
        {
            b +=0;
            b_sig +=0;
            //std::cout << "d = " << std::bitset<64>(d)  << std::endl;

        }
    }
}


void encodekmer_rev(const char *s, uint64_t &b, uint64_t &b_sig){

    using namespace std;

    for (int i = 0; i <int(SL); i++)
    {
        b = b << 1;
        b_sig = b_sig << 2;

        if (s[i] == 'A'){
            b += 4294967297;
            b_sig += 3;
            //std::cout << "d = " << std::bitset<64>(d)  << std::endl;
        }

        else if (s[i] == 'C'){
            b +=4294967296;
            b_sig +=2;
            //std::cout << "d = " << std::bitset<64>(d)  << std::endl;
        }

        else if (s[i] == 'G'){
            b +=1;
            b_sig +=1;
            //std::cout << "d = " << std::bitset<64>(d)  << std::endl;
        }

        else
        {
            b +=0;
            b_sig +=0;
            //std::cout << "d = " << std::bitset<64>(d)  << std::endl;

        }
    }
}



void update_kmer(const char *s, uint64_t &b, uint64_t &b_sig){
    using namespace std;

    b = b << 1;
    b_sig = b_sig << 2;

    //create a mask that has set 32 bit only and set bit 32 to 0
    uint64_t mask = 4294967297;
    b = b & ~mask;

    if (s[0] == 'T')
    {
        b += 4294967297;
        b_sig += 3;

    } else if (s[0] == 'G')
    {
        b += 4294967296;
        b_sig += 2;
    }
    else if (s[0] == 'C')
    {
        b += 1;
        b_sig += 1;
    }
    else
    {
        b += 0;
        b_sig += 0;
    }

}


void update_kmer_rev(const char *s, uint64_t &b, uint64_t &b_sig){
    using namespace std;

    b = b << 1;
    b_sig = b_sig << 2;

    //create a mask that has set 32 bit only and set bit 32 to 0
    uint64_t mask = 4294967297;
    b = b & ~mask;

    if (s[0] == 'A')
    {
        b += 4294967297;
        b_sig += 3;

    } else if (s[0] == 'C')
    {
        b += 4294967296;
        b_sig += 2;
    }
    else if (s[0] == 'G')
    {
        b += 1;
        b_sig += 1;
    }
    else
    {
        b += 0;
        b_sig += 0;
    }

}



bool hd(uint64_t x, uint64_t y, int m){

    //cout << "here" << endl;

    uint64_t z = x ^ y;
    //uint32_t z1 = z;
    uint32_t z2 = z >> 32;
    uint32_t z1 = z;
    uint32_t zc = z1 | z2;

    int ans = 0;
    ans = __builtin_popcount(zc);

    //cout << ans << endl;

    if (ans > m){
        return false;
    }


    return true;
}


vector<string> list_dir(const char *path) {
    vector<string> userString;
    struct dirent *entry;
    DIR *dir = opendir(path);

    while ((entry = readdir(dir)) != NULL) {
        string temp = string(path) + "/" + entry->d_name;
        userString.push_back(temp);
    }
    closedir(dir);

    return(userString);
}

uint64_t encodekmer_bits(uint64_t val, vector<int8_t> shifts, vector<int8_t> bits_to_grab)
{

    using namespace std;

    uint64_t res = 0;
    int i = 0;

    while ( shifts[i]!=-1)
    {
        val = val<<shifts[i];
        asm("shld %b3, %2, %0": "=rm" (res): "0" (res), "r" (val), "ic" ( bits_to_grab[i]) : "cc");

        i++;
    }

    //std::cout << "res, " << std::bitset<32>(res) << "\n";

    return uint64_t(res);
}



uint64_t encodekmer_bits_rev(const char *s, vector<int> pos){
    using namespace std;
    uint64_t d = 0;

    for (int i = 0; i < pos.size(); i++)
    {
        d = d << 2;

        if (s[pos[i]] == 'A'){
            d +=3;
        }

        else if (s[pos[i]] == 'C'){
            d +=2;
        }

        else if (s[pos[i]] == 'G'){
            d +=1;
        }

        else{
            d +=0;
        }

    }

    return d;
}

uint8_t get_encid(uint64_t sind, uint8_t enc_arr_id[])
{

    //cout << sind << endl;
    uint64_t eind = sind >> 3;
    //cout << eind << endl;
    uint64_t ebit = sind % 8;
    //cout << ebit << endl;

    uint8_t enc_arr_indict;
    enc_arr_indict = (enc_arr_id[eind] >> (7-ebit))&1;

    //cout << unsigned(enc_arr_indict) << endl;


    return(enc_arr_indict);
}

