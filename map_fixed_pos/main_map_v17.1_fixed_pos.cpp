// main.cpp -- defining, prototyping, and calling a function

// internal version: main_map_v17.1.cpp

// features:
// change signatures from 64 to 32 bit int
// replace map with arrays
// dump entire array into file without for loop
// get tag from front of the kmer
// based on version 14.1 plus adding two encoding arrays
// collapse tag array indicators into single uint64_t value

// to compile:
// g++-9 main_map_v15.1.cpp -std=c++11 -O3 -o main_map


// to run
// ./main_map -i G000016385_neighbors_copy2.fa -d 3  -o G000016385_nbr_map
// ./main_map -i G000016385_minimized -d 3  -o G000016385_nbr_map
// ./main_map -i G000399765_fna_neighborhood_k32C.fa -d 3 -o G000399765_nbr_map
// ./main_map -i G000307305_dump.fa -d 3 -o G000307305_nbr_map
// ./main_map -i k32C_af_mininimization.fa -o G000307305_nbr_map

// can be run without -d argument

// works with hamming_search 6.2


#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iterator>
#include <vector>
#include <random>
#include <map>
#include <set>
#include <unordered_set>
#include <string>
#include <cstring>
#include <algorithm>
#include <bitset>
#include <unordered_map>
#include <time.h>
#include <utility>
#include <chrono>
#include <new>
#include <numeric>
#include <sys/stat.h>
#include <sys/types.h>
#include <tuple>

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>


#define SL 32
#define SIGS_COLMN 7

#define SIGF_CHUNKS 24
#define TAGF_CHUNKS 24
#define ENCF_CHUNKS 24


using namespace std;

// prototypes

void encodekmer(const char s[], uint64_t &b, uint64_t &b_sig);
void encodekmer_rev(const char *s, uint64_t &b, uint64_t &b_sig);
uint64_t encodekmer_bits(uint64_t val, vector<int8_t> shifts, vector<int8_t> bits_to_grab);
bool hd(uint64_t x, uint64_t y, int m);

uint64_t file_read( istream & is, vector <char> & buff );
uint64_t count_lines( const vector <char> & buff, int sz );



int main(int argc, char *argv[]) {



    //typedef unsigned long long int uint64_t;

    using namespace std;


    auto start = chrono::steady_clock::now();
    srand(time(NULL));


    char *ivalue = NULL;
    uint64_t dvalue = 3;
    char *ovalue = NULL;
    //int index;
    int c;

    opterr = 0;

    while ((c = getopt (argc, argv, "i:d:o:")) != -1)
        switch (c)
        {
            case 'i':
                ivalue = optarg;
                //printf("Option -i has arg: %s\n", optarg);
                break;
            case 'd':
                dvalue = atoi(optarg); // at the moment p=3 & is fixed so dvalue is not being used
                break;
                //printf("Option -d has arg: %s\n", optarg);

            case 'o':
                ovalue = optarg;
                //printf("Option -o has arg: %s\n", optarg);
                break;
            case '?':
                if (optopt == 'i')
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);

                else if (optopt == 'o')
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


    printf ("-- Arguments supplied are \nFilename %s \nOutput %s\n", ivalue, ovalue);

//    for (index = optind; index < argc; index++)
//        printf ("Non-option argument %s\n", argv[index]);



//    cout << ivalue << endl;
//    cout << dvalue << endl;
//    cout << ovalue << endl;



    uint64_t kmer_count = 0;

    if (argc <= 7)
    {
       // printf("-- Arguments supplied are \n");

        //cout << "buffer\n";
        const int SZ = 1024 * 1024;
        vector <char> buff( SZ );

        // open input file
        string input_fin = ivalue;
        ifstream ifs(input_fin);

        if (! ifs.is_open()) {
            // show message:
            std::cout << "Error opening file \n";
            exit(0);
        }

        while( uint64_t cc = file_read( ifs, buff ) )
        {
            kmer_count += count_lines( buff, cc );
        }

        ifs.close();

    }
    else if (argc > 7)
    {
        printf("Too many arguments supplied.\n");
        exit(0);
    }


    string input_fin = ivalue;
    ifstream fin(input_fin);
    if (! fin.is_open()) {
        // show message:
        std::cout << "Error opening file \n";
        exit(0);
    }

    uint64_t p = 3;

//    if ((dvalue <= 0) || (dvalue > 3))
//    {
//        p = 3;
//    }
//    else
//    {
//        //p = dvalue;
//    }

    //int SL = 32;
    uint64_t L = 10;
    float alpha = 0.95;
    uint64_t K = round(log(1 - (pow((1 - alpha), (1 / float(L))))) / log(1 - float(p) / float(SL)));
    K = 15;

    //cout << "Filename " << input_fin << endl;
    cout << "k-mer count = " << kmer_count << endl;
    cout << "p = " << p << '\n';
    cout << "SL = " << int(SL) << '\n';

    L = 2;   // used for local testing
    cout << "L = " << L << '\n';
    cout << "alpha = " << alpha << '\n';
    cout << "Using K = " << K << '\n';





    uint64_t tag_size = 2;  // tag size in bits
    uint64_t partitions = (pow (2, tag_size));;
    uint64_t sigs_row_count = (pow (2, (2*K)-tag_size));
    //cout << sigs_row_count << endl;


    // allocate signature array
    uint32_t * sigs_arr;
    uint64_t sigs_arr_size = sigs_row_count*SIGS_COLMN*partitions*L;

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
    uint8_t * sigs_indicator_arr;
    uint64_t sigs_indicator_arr_size = sigs_row_count*L;

    try
    {
        sigs_indicator_arr = new uint8_t [sigs_indicator_arr_size];
        cout << "-- Done indicator allocation " << endl;
    }
    catch (std::bad_alloc& ba)
    {
        std::cerr << "bad_alloc caught: " << ba.what() << '\n';
    }


    // initialize indicator array to 0
    for(uint64_t enc=0; enc< sigs_indicator_arr_size; enc++)
    {
        sigs_indicator_arr[enc]= 0;

    }



    // allocate tag array
    int8_t * tag_arr;
    uint64_t tag_arr_size = sigs_row_count*SIGS_COLMN*partitions*L;

    try
    {
        tag_arr = new int8_t [tag_arr_size];
        cout << "-- Done tag allocation " << endl;
    }
    catch (std::bad_alloc& ba)
    {
        std::cerr << "bad_alloc caught: " << ba.what() << '\n';
    }


    // initialize tag array to -1

//    for(uint64_t j = 0; j < sigs_row_count*SIGS_COLMN*partitions*L; j++)
//    {
//        tag_arr[j] = -1;
//    }

    uint64_t MAX_ENC_CNT = 4000000000;


    // allocate encoding array
    uint64_t * encode_arr; // AT
    uint64_t * encode_arr_1; // CG
    uint64_t encode_arr_size;

    if (kmer_count < MAX_ENC_CNT)
    {
        encode_arr_size = kmer_count;
    }
    else
    {
        encode_arr_size = MAX_ENC_CNT+5;
    }



    try
    {
        encode_arr = new uint64_t [encode_arr_size];
        cout << "-- Done encoding array 0 allocation " << endl;
    }
    catch (std::bad_alloc& ba)
    {
        std::cerr << "bad_alloc caught: " << ba.what() << '\n';
    }

    try
    {
        encode_arr_1 = new uint64_t [encode_arr_size];
        cout << "-- Done encoding array 1 allocation " << endl;
    }
    catch (std::bad_alloc& ba)
    {
        std::cerr << "bad_alloc caught: " << ba.what() << '\n';
    }


    // allocate enc array id array
    uint8_t * enc_arr_id;
    uint64_t enc_arr_id_size = sigs_arr_size;


    try
    {
        enc_arr_id = new uint8_t [enc_arr_id_size];
        cout << "-- Done enc id array allocation " << endl;
    }
    catch (std::bad_alloc& ba)
    {
        std::cerr << "bad_alloc caught: " << ba.what() << '\n';
    }

    //initialize enc_arr_id to 0
    for(uint64_t j = 0; j < enc_arr_id_size; j++)
    {
        enc_arr_id[j] = 0;
    }



    // write map to file
    string map_name = ovalue;

    // create a directory
    // linux dependant function
    const int dir_err = mkdir(map_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    if (-1 == dir_err)
    {
        //printf("Error creating directory!n");
        cerr << "Error :  " << strerror(errno) << endl;
        exit(1);
    }


    // generate mask
    uint64_t l, k, n;


    uint64_t included_kmers_counter[L];

    vector <vector<int> > positions (L, vector<int> (K));

    srand(time(NULL));

    for (l = 0; l < L; l++) {
        vector<int> rand_num;

        // initialize lost kmers array
        included_kmers_counter[l] = 0;


        for (k = 0; k < K; k++) {
            n = rand() % int(SL);

            if (count(rand_num.begin(), rand_num.end(), n)) {
                k -= 1;
            }
            else {
                rand_num.push_back(n);
            }

        }

        cout << "Positions for l : " << l << " >> ";

        sort(rand_num.begin(), rand_num.end(), std::greater<int>());
        for (int j =0; j < K; j++)
        {
            positions[l][j] = rand_num[j];
            cout << positions[l][j] << ", ";
        }

        cout << endl;

        rand_num.clear();
    }

    // mask from 15.2 with k=14 experiment
    //positions [0] = {31, 29, 28, 26, 25, 22, 18, 16, 10, 8, 6, 4, 2, 0};
    //positions [1] = {30, 28, 27, 26, 21, 19, 18, 17, 16, 13, 8, 4, 3, 0};

    // mask from 15.2 with k=15 experiment
    positions [0] = {30, 28, 25, 24, 22, 21, 20, 17, 14, 13, 11, 5, 3, 2, 1};
    positions [1] = {30, 29, 26, 24, 20, 19, 16, 14, 13, 12, 10, 7, 3, 1, 0};
    //positions [2] = {30, 29, 26, 23, 21, 20, 19, 18, 15, 11, 9, 7, 5, 3, 2};

    for (l = 0; l < L; l++)
    {
        cout << "New positions for l : " << l << " >> ";

        for (int j =0; j < K; j++)
        {
            cout << positions[l][j] << ", ";
        }

        cout << endl;
    }



    vector<vector<int8_t>> shifts;
    vector<vector<int8_t>> grab_bits;


    for (l = 0; l < L; l++)
    {
        // construct a vector of int
        vector<int8_t> v;
        vector<int8_t> g;
        int lp = 31;
        int jp = 0;


        for (int p = 0; p < K; p++)
        {

            if (p == 0)
            {
                v.push_back((lp-positions[l][p])*2);
                lp = positions[l][p];
                jp +=2;

            }
            else if (positions[l][p-1] - positions[l][p] !=1)
            {
                v.push_back((lp - positions[l][p])*2);
                lp = positions[l][p];
                g.push_back(jp);
                jp = 2;
            }
            else
            {
                jp+=2;
            }

        }

        g.push_back(jp);

        // push back above one-dimensional vector
        v.push_back(-1);
        g.push_back(-1);

        shifts.push_back(v);
        grab_bits.push_back(g);
    }



//    int q  = 0;
//    while (shifts[0][q] !=-1)
//    {
//        cout << "shifts " << unsigned(shifts[0][q]) << endl;
//        cout << "grab_bits " << unsigned(grab_bits[0][q]) << endl;
//        q++;
//    }


    string line;
    uint64_t li = 0;
    uint64_t encli_0 = 0;
    uint64_t encli_1 = 0;

    uint64_t b;
    uint64_t b_sig;
    bool kmerwritten;

    uint64_t sig_hash;
    int8_t tag;
    uint64_t big_sig_hash;

    uint8_t enc_array_indict = 0;

    // tag size is number of bits extracted from the right of the signature value
    //masks if tag = 2
//    uint8_t  tag_mask = 3;               // 3
//    uint32_t big_sig_mask = 67108863;    // 3FFFFFF

    // compute mask values
    uint64_t  tag_mask = 0;

    for (int i = 0; i < tag_size; i ++)
    {
        tag_mask += (pow (2, i)); // 111 if tag = 3
    }
    //cout << "tag mask " << unsigned(tag_mask) << endl;


    uint64_t big_sig_mask = 0;

    for (int i = 0; i < ((2*K)-tag_size); i++)
    {
        big_sig_mask += (pow (2, i));
    }
    //cout << "big sig mask " << big_sig_mask << endl;




    while (std::getline(fin, line)) {

        if (line.rfind('>', 0) != 0) {
            const char *cline = line.c_str();

            b = 0;
            b_sig = 0;
            kmerwritten = false;


            encodekmer(cline, b, b_sig);

//            std::cout << "b = " << std::bitset<64>(b)  << std::endl;
//            std::cout << "b_sig = " << std::bitset<64>(b_sig)  << std::endl;



            if (li % 1000000 == 0) {
                cout << "-- Encoding  " << li << endl;
            }

            // set encoding array id
            enc_array_indict = rand()%2;
            //enc_array_indict = 1;


            // check if maximum number of encodings (0xFFFFFFFF or 4294967295) is not exceeded
            if ((encli_0 >= MAX_ENC_CNT) || (encli_1 >= MAX_ENC_CNT))
            {
                break;
            }

                for (l = 0; l < L; l++) {

                   sig_hash = encodekmer_bits(b_sig, shifts[l], grab_bits[l]);
                   //std::cout << "sig_hash = " << std::bitset<32>(sig_hash)  << std::endl;

                   // get first 2 bits of signature (effectively bits 28 - 27) of 32 bit encoding as tag
                   tag = (sig_hash >>((2*K)-tag_size))& tag_mask;
                   //std::cout << "tag = " << std::bitset<8>(tag)  << std::endl;


                   // get last 26 bits of signature (effectively bits 26 -1 ) of 32 bit encoding as sigs row number
                   big_sig_hash = sig_hash& big_sig_mask;
                   //std::cout << "big_sig = " << std::bitset<32>(big_sig_hash)  << std::endl;

                   //cout << sig_hash << endl;
                   //cout << "tag " << unsigned(tag) << endl;
                   //cout << big_sig_hash << endl;


                   // check is row space is available for forward kmer
                   if (sigs_indicator_arr[sigs_row_count*l + big_sig_hash]< ((uint64_t)SIGS_COLMN)*partitions)
                   {
                        // populate tag array
                        tag_arr[sigs_row_count*((uint64_t)SIGS_COLMN)*partitions*l + big_sig_hash*((uint64_t)SIGS_COLMN)*partitions +
                                sigs_indicator_arr[sigs_row_count*l + big_sig_hash]] = tag;

                        if (enc_array_indict == 0 )  // first letter is either A or C
                        {
                            // populate sigs array
                            sigs_arr[sigs_row_count*((uint64_t)SIGS_COLMN)*partitions*l + big_sig_hash*((uint64_t)SIGS_COLMN)*partitions +
                                     sigs_indicator_arr[sigs_row_count*l + big_sig_hash]]= encli_0;

                        }
                        else if (enc_array_indict == 1 ) // first letter is either G or T
                        {
                            // populate sigs array
                            sigs_arr[sigs_row_count*((uint64_t)SIGS_COLMN)*partitions*l + big_sig_hash*((uint64_t)SIGS_COLMN)*partitions +
                                     sigs_indicator_arr[sigs_row_count*l + big_sig_hash]]= encli_1;
                        }

                        // populate enc array id array
                        enc_arr_id[sigs_row_count*((uint64_t)SIGS_COLMN)*partitions*l + big_sig_hash*((uint64_t)SIGS_COLMN)*partitions +
                                   sigs_indicator_arr[sigs_row_count*l + big_sig_hash]]= enc_array_indict;

                        // increment indicator array
                        sigs_indicator_arr[sigs_row_count*l + big_sig_hash] +=1;

                        // increment kmer counter
                        included_kmers_counter[l]+=1;

                        kmerwritten = true;
                    }


                    //cout <<  *(sigs_arr + l*kmer_count + li) << endl;
                }



            if (kmerwritten)
            {
                if  (enc_array_indict == 0)
                {
                    encode_arr[encli_0] = b;
                    encli_0 +=1;
                }
                else
                {
                    encode_arr_1[encli_1] = b;
                    encli_1 +=1;
                }
            }

            li += 1;

        }
    }

    // update kmer count
    kmer_count = encli_0 + encli_1;

    fin.close();


    // recording end time.
    auto end = chrono::steady_clock::now();
    cout << "-- Done hashing.  Now writing. Time so far: "
         << chrono::duration_cast<chrono::seconds>(end - start).count() << " seconds" << endl;




    // allocate new tag array
    int8_t * new_tag_arr;

    uint64_t new_tag_arr_size = sigs_row_count*partitions*L;

    try
    {
        new_tag_arr = new int8_t [new_tag_arr_size];
        cout << "-- Done new tag allocation " << endl;
    }
    catch (std::bad_alloc& ba)
    {
        std::cerr << "bad_alloc caught: " << ba.what() << '\n';
    }

    // initialize new tag array to -1
    for(uint64_t j = 0; j < new_tag_arr_size; j++){
        new_tag_arr[j] = -1;
    }



    // allocate new enc_ind_array
    uint8_t * new_encid_arr;
    //uint64_t new_encid_arr_clmn = SIGS_COLMN*partitions/8;
    uint64_t increment = 8;
    uint64_t new_encid_arr_sz = ceil(sigs_arr_size/increment);

    try
    {
        new_encid_arr = new uint8_t [new_encid_arr_sz];
        cout << "-- Done new enc id allocation " << endl;
    }
    catch (std::bad_alloc& ba)
    {
        std::cerr << "bad_alloc caught: " << ba.what() << '\n';
    }

//    // initialize new encid array to 0
//    for(uint64_t j = 0; j < new_encid_arr_sz; j++){
//        new_encid_arr[j] = 0;
//    }





    // sort sigs and tag arrays
    // sort uint32_t array b[] according to the order defined by a[]
    for(l = 0; l < L; l++)
    {

        for (uint64_t j =0; j < sigs_row_count; j++)
        {

            if (sigs_indicator_arr[sigs_row_count*l + j]>0)
            {

                vector<int8_t> tag_clmn_count_vec (partitions, 0);

                uint8_t n = sigs_indicator_arr[sigs_row_count*l + j];
                //cout << unsigned(n) << endl;



                //tuple<int8_t, uint32_t, uint8_t> pairt[n];
                //triple pairt[n];

                //vector <tuple<int8_t, uint32_t, uint8_t> > pairt (n);

                vector <tuple<int8_t, uint32_t, uint8_t> > pairt;



                // Storing the respective array elements in pairs.
                for (uint64_t i = 0; i < n; i++)
                {
                    //cout << i << endl;
//                    pairt[i].tag_val = tag_arr[sigs_row_count*SIGS_COLMN*partitions*l + j*SIGS_COLMN*partitions+i];
//                    pairt[i].sig_val = sigs_arr[sigs_row_count*SIGS_COLMN*partitions*l + j*SIGS_COLMN*partitions+i];
//                    pairt[i].encoding_id = enc_arr_id[sigs_row_count*SIGS_COLMN*partitions*l + j*SIGS_COLMN*partitions+i];
//                    //std::cout << "pair.second = " << std::bitset<64>(pairt[i].second)  << std::endl;
//
//                    tag_clmn_count_vec[pairt[i].tag_val] +=1;

                    int8_t tag_val = tag_arr[sigs_row_count*((uint64_t)SIGS_COLMN)*partitions*l + j*((uint64_t)SIGS_COLMN)*partitions+i];
                    uint32_t sig_val = sigs_arr[sigs_row_count*((uint64_t)SIGS_COLMN)*partitions*l + j*((uint64_t)SIGS_COLMN)*partitions+i];
                    uint8_t encoding_id = enc_arr_id[sigs_row_count*((uint64_t)SIGS_COLMN)*partitions*l + j*((uint64_t)SIGS_COLMN)*partitions+i];

                    //pairt[n] = (make_tuple(tag_val, sig_val, encoding_id));
                    pairt.push_back(make_tuple(tag_val, sig_val, encoding_id));


                    tag_clmn_count_vec[tag_val] +=1;

                    //tag_arr[sigs_row_count*SIGS_COLMN*partitions*l + j*SIGS_COLMN*partitions+i] = -1;

                }

//                if ( n == 24) {
//                    for (int k = 0; k < 24; k++) {
//                        cout << unsigned(pairt[k].tag_val) << " : " << unsigned(pairt[k].encoding_id) << ", ";
//                    }
//                    cout << endl;
//                }

                // Sorting the pair array.
                //sort(pairt.begin(), pairt.end());

                //sort(pairt, pairt + n, compare);
                sort(pairt.begin(), pairt.end());

                //std::sort(pairt, pairt + n, [] (auto t1, auto t2) {return t1.tag_val < t2.tag_val;});

//                if ( n == 24)
//                {
//                    for (int k=0; k < 24; k++)
//                    {
//                        cout << unsigned(pairt[k].tag_val) << " : " << unsigned(pairt[k].encoding_id) << ", " ;
//                    }
//
//                    cout << endl;
//                }


                // update new_tag_array
                for (uint64_t r = 0; r < partitions; r++)
                {
                    if (r == 0)
                    {
                        new_tag_arr[sigs_row_count*partitions*l + j*partitions + r] = tag_clmn_count_vec[r];
                        //cout << unsigned(new_tag_arr[sigs_row_count*partitions*l + j*partitions + r]) << endl;
                    }
                    else
                    {
                        new_tag_arr[sigs_row_count*partitions*l + j*partitions + r] =
                                new_tag_arr[sigs_row_count*partitions*l + j*partitions + r-1] + tag_clmn_count_vec[r];
                        //cout << unsigned(new_tag_arr[sigs_row_count*partitions*l + j*partitions + r]) << endl;

                    }


                }

                // Modifying original sig array
                for (uint64_t i = 0; i < n; i++)
                {
//                    sigs_arr[sigs_row_count*SIGS_COLMN*partitions*l + j*SIGS_COLMN*partitions+i] = pairt[i].sig_val;
//                    enc_arr_id[sigs_row_count*SIGS_COLMN*partitions*l + j*SIGS_COLMN*partitions+i] = pairt[i].encoding_id;
                    sigs_arr[sigs_row_count*((uint64_t)SIGS_COLMN)*partitions*l + j*((uint64_t)SIGS_COLMN)*partitions+i] = get<1>(pairt[i]);
                    enc_arr_id[sigs_row_count*((uint64_t)SIGS_COLMN)*partitions*l + j*((uint64_t)SIGS_COLMN)*partitions+i] = get<2>(pairt[i]);
                }

                //fill(tag_clmn_count_vec.begin(), tag_clmn_count_vec.end(), 0);
                //pairt.clear();
                // merge every set of 8 bis into single value in enc_arr_id array


            }

        }

    }


    // populate new encid array
    uint64_t j = 0;

    for (uint64_t i = 0; i < enc_arr_id_size; i += increment)
    {
        uint8_t combo_id = 0;

        for (uint64_t x = 0; x < increment; x++){
            combo_id = combo_id << 1;
            combo_id = combo_id | enc_arr_id[i + x];
        }

        new_encid_arr[j] = combo_id;
        j +=1;

    }




//    if (mkdir(map_name.c_str(), 0777) == -1)
//        cerr << "Error :  " << strerror(errno) << endl;
//
//    else
//        cout << "Directory created";





    // open metadata file

    string map_meta = map_name + "_meta";
    string path = map_name + "/" + map_meta;

    FILE *wfmeta;
    wfmeta = fopen(path.c_str(), "wb");
    if (! wfmeta) {
        cout << "Cannot open file!" << endl;
        return 1;
    }




    // write metadata file

    // write p, L, alpha, K values
    fwrite(&p, sizeof(uint64_t), 1, wfmeta);
    //int SL_val = SL;
    //wf.write((char *)&SL_val, sizeof(int8_t));
    fwrite(&L, sizeof(uint64_t), 1, wfmeta);
    fwrite(&alpha, sizeof(float), 1, wfmeta);
    fwrite(&K, sizeof(uint64_t), 1, wfmeta);


    fwrite(&sigs_arr_size, sizeof(uint64_t), 1, wfmeta);
    fwrite(&new_tag_arr_size, sizeof(uint64_t), 1, wfmeta);
    //cout << "kmer count " << kmer_count << endl;

    // write k-mer count
    fwrite(&kmer_count, sizeof(uint64_t), 1, wfmeta);
    cout << "kmer count " << kmer_count << endl;
    fwrite(&encli_0, sizeof(uint64_t), 1, wfmeta);
    cout << "kmer count array 0 " << encli_0 << endl;
    fwrite(&encli_1, sizeof(uint64_t), 1, wfmeta);
    cout << "kmer count array 1 " << encli_1 << endl;

    fwrite(&new_encid_arr_sz, sizeof(uint64_t), 1, wfmeta);



    // write mask array and output real count k-mers included in DB

    for (l=0; l < L; l++)
    {
        cout << "k-mers included l " << l <<  " : " << included_kmers_counter[l] << endl;


        int size = shifts[l].size();
        //cout << "size " << size << endl;
        fwrite(&size, sizeof(int8_t), 1, wfmeta);

        for (int s =0; s < size; s++){
            fwrite(&shifts[l][s], sizeof(int8_t), 1, wfmeta);
            //cout << signed(shifts[l][s]) << endl;
        }



        for (int s = 0; s < size; s++){
            fwrite(&grab_bits[l][s], sizeof(int8_t), 1, wfmeta);
            //cout << signed(grab_bits[l][s]) << endl;
        }

    }


    // write sig columns counts and  sig row count
    uint64_t sig_clmn_cnt = SIGS_COLMN;
    fwrite(&sig_clmn_cnt, sizeof(uint64_t), 1, wfmeta);
    fwrite(&sigs_row_count, sizeof(uint64_t), 1, wfmeta);


    // write tag size, partition count, tag mask and big_sig_mask
    fwrite(&tag_size, sizeof(uint64_t), 1, wfmeta);
    fwrite(&partitions, sizeof(uint64_t), 1, wfmeta);
    fwrite(&tag_mask, sizeof(uint64_t), 1, wfmeta);
    fwrite(&big_sig_mask, sizeof(uint64_t), 1, wfmeta);


    // write file chunk information
    uint8_t sigf_chunks = SIGF_CHUNKS;
    uint8_t tagf_chunks = TAGF_CHUNKS;
    uint8_t encf_chunks = ENCF_CHUNKS;
    fwrite(&sigf_chunks, sizeof(uint8_t), 1, wfmeta);
    fwrite(&tagf_chunks, sizeof(uint8_t), 1, wfmeta);
    fwrite(&encf_chunks, sizeof(uint8_t), 1, wfmeta);



    // initialize counter variables for sigs, indicator and tags arrays to output counts
    uint64_t total_sigs_written = 0;
    //uint64_t total_indicators_written = 0;
    uint64_t total_tags_written = 0;
    uint64_t total_members_written = 0;
    uint64_t total_members_written_1 = 0;
    uint64_t total_encid_written = 0;



    // compute members in file chunks
    vector<uint64_t> sig_chunk_cnts(SIGF_CHUNKS, 0);
    vector<uint64_t> tag_chunk_cnts(TAGF_CHUNKS, 0);
    vector<uint64_t> enc_chunk_cnts(ENCF_CHUNKS, 0);
    vector<uint64_t> enc_chunk_cnts_1(ENCF_CHUNKS, 0);
    vector<uint64_t> encid_chunk_cnts(SIGF_CHUNKS, 0);


    for (int m =0; m < SIGF_CHUNKS-1; m++){
        sig_chunk_cnts[m] = round(sigs_arr_size/SIGF_CHUNKS);
        total_sigs_written += sig_chunk_cnts[m];
    }
    sig_chunk_cnts[SIGF_CHUNKS-1] = sigs_arr_size - total_sigs_written;



    for (int m =0; m < TAGF_CHUNKS-1; m++){
        tag_chunk_cnts[m] = round(new_tag_arr_size/TAGF_CHUNKS);
        total_tags_written += tag_chunk_cnts[m];
    }
    tag_chunk_cnts[TAGF_CHUNKS-1] = new_tag_arr_size - total_tags_written;



    for (int m =0; m < ENCF_CHUNKS-1; m++){
        enc_chunk_cnts[m] = round(encli_0/ENCF_CHUNKS);
        total_members_written += enc_chunk_cnts[m];
    }
    enc_chunk_cnts[ENCF_CHUNKS-1] = encli_0 - total_members_written;


    for (int m =0; m < ENCF_CHUNKS-1; m++){
        enc_chunk_cnts_1[m] = round(encli_1/ENCF_CHUNKS);
        total_members_written_1 += enc_chunk_cnts_1[m];
    }
    enc_chunk_cnts_1[ENCF_CHUNKS-1] = encli_1 - total_members_written_1;


    for (int m =0; m < SIGF_CHUNKS-1; m++){
        encid_chunk_cnts[m] = round(new_encid_arr_sz/SIGF_CHUNKS);
        total_encid_written += encid_chunk_cnts[m];
    }
    encid_chunk_cnts[SIGF_CHUNKS-1] = new_encid_arr_sz - total_encid_written;



    total_sigs_written = 0;
    total_tags_written = 0;
    total_members_written = 0;
    total_members_written_1 = 0;
    total_encid_written = 0;

    // check total counts and write to file
    for (int m = 0; m < SIGF_CHUNKS; m ++)
    {
        fwrite(&sig_chunk_cnts[m], sizeof(uint64_t), 1, wfmeta);

        // open sig file
        string map_sig = map_name + "_sig" + to_string(m);
        path = map_name + "/" + map_sig;

        FILE *wf;
        wf = fopen(path.c_str(), "wb");

        if (!wf) {
            cout << "Cannot open file!" << endl;
            return 1;
        }

            // write sigs
            total_sigs_written += fwrite(sigs_arr+total_sigs_written, sizeof(uint32_t), sig_chunk_cnts[m], wf);
            fclose(wf);
        }






    // check total counts and write to file
    for (int m = 0; m < TAGF_CHUNKS; m++)
    {
        fwrite(&tag_chunk_cnts[m], sizeof(uint64_t), 1, wfmeta);

        // open tag file
        string map_tag = map_name + "_tag" + to_string(m);
        path = map_name + "/" + map_tag;


        FILE *wftag;
        wftag = fopen(path.c_str(), "wb");
        if (! wftag) {
            cout << "Cannot open file!" << endl;
            return 1;
        }

        // write tags
        total_tags_written += fwrite(new_tag_arr + total_tags_written, sizeof(int8_t), tag_chunk_cnts[m], wftag);
        fclose(wftag);
    }




    // check total counts and write to file
    for (int m = 0; m < ENCF_CHUNKS; m ++)
    {
        fwrite(&enc_chunk_cnts[m], sizeof(uint64_t), 1, wfmeta);

        // open encoding file

        string map_enc = map_name + "_enc" + to_string(m);
        path = map_name + "/" + map_enc;


        FILE *wfenc;
        wfenc = fopen(path.c_str(), "wb");
        if (! wfenc) {
            cout << "Cannot open file!" << endl;
            return 1;
        }

        // write encoding array
        total_members_written += fwrite(encode_arr + total_members_written, sizeof(uint64_t), enc_chunk_cnts[m], wfenc);
        fclose(wfenc);
    }


    // check total counts and write to file
    for (int m = 0; m < ENCF_CHUNKS; m ++)
    {
        fwrite(&enc_chunk_cnts_1[m], sizeof(uint64_t), 1, wfmeta);

        // open encoding file

        string map_ence = map_name + "_ence" + to_string(m);
        path = map_name + "/" + map_ence;


        FILE *wfenc;
        wfenc = fopen(path.c_str(), "wb");
        if (! wfenc) {
            cout << "Cannot open file!" << endl;
            return 1;
        }

        // write encoding array
        total_members_written_1 += fwrite(encode_arr_1 + total_members_written_1, sizeof(uint64_t), enc_chunk_cnts_1[m], wfenc);
        fclose(wfenc);
    }


    // check total encid and write to file
    for (int m = 0; m < SIGF_CHUNKS; m ++)
    {
        fwrite(&encid_chunk_cnts[m], sizeof(uint64_t), 1, wfmeta);

        // open sig file
        string map_encid = map_name + "_encid" + to_string(m);
        path = map_name + "/" + map_encid;

        FILE *wf;
        wf = fopen(path.c_str(), "wb");

        if (!wf) {
            cout << "Cannot open file!" << endl;
            return 1;
        }

        // write sigs
        total_encid_written += fwrite(new_encid_arr+total_encid_written, sizeof(uint8_t), encid_chunk_cnts[m], wf);
        fclose(wf);
    }

    //cout << unsigned(SIGS_COLMN) << " " << sigs_row_count << " " << unsigned(tag_size) << " " << unsigned(partitions) << " " << unsigned(tag_mask) << " " << big_sig_mask << endl;
    //total_indicators_written = fwrite(sigs_indicator_arr, sizeof(uint8_t), sigs_row_count*L, wf);


    cout << "Sigs written  " <<  ": " << total_sigs_written << endl;
    //cout << "Indicators written  " << ": " << total_indicators_written << endl;
    cout << "Tags written  " << ": " << total_tags_written << endl;
    cout << "Encodings written " << ": " << total_members_written + total_members_written_1 << endl;
    cout << "Encodings written to array 0 " << ": " << total_members_written << endl;
    cout << "Encodings written to array 1 " << ": " << total_members_written_1 << endl;
    cout << "Enc id written " << ": " << total_encid_written << endl;


    //fclose(wf);
    fclose(wfmeta);
    //fclose(wftag);
    //fclose(wfenc);


    end = chrono::steady_clock::now();
    cout << "-- Done making map.  Time so far: " << chrono::duration_cast<chrono::seconds>(end - start).count()
         << " seconds" << endl;


    // output map information
    cout << "Columns = " << unsigned(SIGS_COLMN) << endl;
    cout << "Partitions = " << unsigned(partitions) << endl;
    cout << "Sigs row count = " << sigs_row_count << '\n';
    cout << "Tag size = " << unsigned(tag_size) << '\n';
    cout << "Tag mask = " << unsigned(tag_mask) << '\n';
    cout << "Big sig mask = " << big_sig_mask << '\n';


    //compute the usage of rows in signature matrix to determine how many positions in each row are used

    // row count vector
    vector<uint64_t> sig_row_count_vec (SIGS_COLMN * partitions+1, 0);


    // traverse indicator array to compute filled positions
    for (int l =0; l < L; l++) {

        for (uint64_t r = 0; r < sigs_row_count; r++) {

            sig_row_count_vec[sigs_indicator_arr[sigs_row_count*l + r]] += 1;

        }
        // output counts for a given l

        //cout << sigs_row_count << endl;

        for (int s = 0; s < SIGS_COLMN * partitions+1; s++)
        {

            cout << "l : " << l << " -- Count of rows with positions filled " << s <<  " : " << sig_row_count_vec[s] << " : ";
            cout << std::fixed << std::setprecision(6) << (double)sig_row_count_vec[s]/sigs_row_count << endl;
        }

        cout << "l : " << l << " -- Total populated rows " << ": " << std::accumulate( sig_row_count_vec.begin(), sig_row_count_vec.end(),
                                                                                       decltype( sig_row_count_vec)::value_type(0)) << endl;
        // reset counts to 0
        fill(sig_row_count_vec.begin(), sig_row_count_vec.end(), 0);

    }



    delete [] sigs_arr ; // remember to delete array when I am done
    delete [] sigs_indicator_arr ; // remember to delete array when I am done
    delete [] tag_arr ; // remember to delete array when I am done
    delete [] encode_arr ; // remember to delete array when I am done
    delete [] encode_arr_1 ; // remember to delete array when I am done
    delete [] enc_arr_id ;
    delete [] new_tag_arr ;



    end = chrono::steady_clock::now();
    cout << "-- Done writing.  Time so far: " << chrono::duration_cast<chrono::seconds>(end - start).count()
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




uint64_t file_read( istream & is, vector <char> & buff ) {
    is.read( &buff[0], buff.size() );
    return is.gcount();
}

uint64_t count_lines( const vector <char> & buff, int sz ) {
    uint64_t newlines = 0;
    const char * p = &buff[0];
    for ( int i = 0; i < sz; i++ ) {
        if ( p[i] == '>' ) {
            newlines++;
        }
    }
    return newlines;
}

