// main.cpp -- defining, prototyping, and calling a function

// features:
// change signatures from 64 to 32 bit int
// replace map with arrays
// dump entire array into file without for loop

// to compile:
// g++-9 minimization_v3.0.cpp -std=c++11 -o main_minimization

// to run:
// ./main_minimization -i G000016385_neighbors_copy2.fa -o G000016385_minimized
// ./main_minimization -i G000399765_fna_neighborhood_k32C.fa -o G000399765_minimized
// ./main_minimization -i G000307305_dump.fa -o G000307305_minimized
// ./main_minimization -i k35C_bef_mininimization.fa -o k32C_af_mininimization.fa

// works with hamming_search 6.2

// need to run jellyfish on the output to collapse kmer minimizers //

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
#include <unordered_set>
#include <limits>


#define ver_num 3.0
#define SL 35
#define MINIMIZER 32
#define SIGS_COLMN 6



using namespace std;

// prototypes

uint64_t encodekmer(const char *s);
string decodekmer(uint64_t d);

uint64_t compute_ascii(const char *s);

uint64_t update_kmer(uint64_t old_b, const char *s);

uint64_t file_read( istream & is, vector <char> & buff );
uint64_t count_lines( const vector <char> & buff, int sz );



int main(int argc, char *argv[]) {

    //typedef unsigned long long int uint64_t;

    using namespace std;

    auto start = chrono::steady_clock::now();
    srand(time(NULL));

    // output version number
    cout << "v." << std::fixed << std::setprecision(1)  << ver_num << endl;



    if (argc <= 5)
    {
        printf("-- Arguments supplied are \n");

    }
    else if (argc > 5)
    {
        printf("Too many arguments supplied.\n");
        exit(0);
    }

    // open input file
    string input_fin = argv[2];
    ifstream fin(input_fin);

    if (! fin.is_open()) {
        // show message:
        std::cout << "Error opening file";
        exit(0);
    }



    // open output file
    ofstream outputFile;
    outputFile.open(argv[4]);

    if (!outputFile.is_open())
    {
        cout << "Cannot open file!" << endl;
        return 1;
    }


    //uint32_t p = 3;
    //int SL = 32;
    //uint32_t L = 10;
    //float alpha = 0.95;
    //uint32_t K = round(log(1 - (pow((1 - alpha), (1 / float(L))))) / log(1 - float(p) / float(SL)));


    cout << "Filename " << input_fin << endl;
    //cout << "k-mer count = " << kmer_count << endl;
    //cout << "p = " << p << '\n';
    cout << "SL = " << int(SL) << '\n';
    cout << "Minimizer length = " << MINIMIZER << '\n';
//    cout << "alpha = " << alpha << '\n';
//    cout << "Using K = " << K << '\n';


    //L = 1;   // used for local testing


    //uint8_t tag_size = 2;  // tag size in bits
    //uint8_t partitions = 4;
    //uint64_t sigs_row_count = (pow (2, (2*K)-tag_size));
    //cout << sigs_row_count << endl;



    string line;
    string name;
    uint64_t kmer_count = 0;
    uint64_t minimizer_count = 0;



    while (std::getline(fin, line)) {

        if  (line.rfind('>', 0) == 0) {
            name = line;
        }

        if (line.rfind('>', 0) != 0) {

            std::string min_str (MINIMIZER, 'Z');


            for (uint64_t i = 0; i < line.length()-MINIMIZER+1; i++) {

                string kmer_str = line.substr(i, int(MINIMIZER));
                //cout << kmer_str << endl;


                if (kmer_str < min_str){
                    min_str = kmer_str;
                }

            }

            if (kmer_count % 1000000 == 0) {
                cout << "-- Minimizing  " << kmer_count << endl;
            }

            // write minimizer for a given kmer to file
            outputFile << name << endl;
            outputFile << min_str << endl;
            minimizer_count +=1;

            kmer_count +=1;




        }
    }


    fin.close();
    outputFile.close();

    //std::cout << "b = " << std::bitset<64>(b)  << std::endl;


    cout << "Total kmers  " <<  ": " << kmer_count << endl;
    cout << "Total minimizers  " << ": " << minimizer_count << endl;
    //cout << "Total members written " << " : " << total_members_written/2 << endl;



    // recording end time
    auto end = chrono::steady_clock::now();
    cout << "-- Done writing.  Time so far: " << chrono::duration_cast<chrono::seconds>(end - start).count()
    << " seconds" << endl;

    return 0;

}



// function definition
uint64_t encodekmer(const char *s){
    using namespace std;
    uint64_t d = 0;

    for (int i = 0; i <int(MINIMIZER); i++)
    {
        d = d << 2;

        if (s[i] == 'T'){
            d +=3;
        }

        else if (s[i] == 'G'){
            d +=2;
        }

        else if (s[i] == 'C'){
            d +=1;
        }

        else{
            d +=0;
        }

//        bitset<32> x(d);
//        cout << x << endl;
    }

    return d;
}


// function definition
string decodekmer(uint64_t d){
    using namespace std;
    string kmer;


    for (int i = 0; i <int(MINIMIZER); i++)
    {
        int val = d & 3;

        if (val == 3){
            kmer += "T";

        }

        else if (val == 2){
            kmer += "G";
        }

        else if (val == 1){
            kmer += "C";
        }

        else{
            kmer += "A";
        }

        d = d >> 2;
//        bitset<32> x(d);
//        cout << x << endl;
    }
    reverse(kmer.begin(), kmer.end());

    return kmer;
}

// function definition
uint64_t compute_ascii(const char *s){
    using namespace std;
    uint64_t d = 0;

    for (int i = 0; i <int(MINIMIZER); i++)
    {
        d += int(s[i]);
    }

    return d;
}

uint64_t update_kmer(uint64_t old_b, const char *s){
    using namespace std;
    //std::cout << "old b = " << std::bitset<64>(old_b )  << std::endl;
    old_b = old_b << 1;


    //create a mask that has set 32 bit only
    uint64_t mask = 4294967297;

    // set bit 32 to 0
    //std::cout << "old b = " << std::bitset<64>(old_b )  << std::endl;
    old_b = old_b & ~mask;
    //std::cout << "old b = " << std::bitset<64>(old_b )  << std::endl;
    //cout << " " << s << endl;


    //std::cout << "old b = " << std::bitset<64>(old_b )  << std::endl;

    //std::cout << "shifted old b = " << std::bitset<64>(old_b )  << std::endl;

    if (s[0] == 'T'){
        old_b +=4294967297;
        //cout << "here " << endl;
        //std::cout << "d = " << std::bitset<64>(d)  << std::endl;
    }

    else if (s[0] == 'G'){
        old_b +=4294967296;
        //std::cout << "d = " << std::bitset<64>(d)  << std::endl;
    }

    else if (s[0] == 'C'){
        old_b +=1;
        //std::cout << "d = " << std::bitset<64>(d)  << std::endl;
    }

    else{
        old_b +=0;
        //std::cout << "d = " << std::bitset<64>(d)  << std::endl;
    }

    //std::cout << "d = " << std::bitset<64>(old_b)  << std::endl;


    return old_b;
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

