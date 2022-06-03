#include "BipartiteGraph.h"
#include "StableMarriage.h"
#include "Convert_HR_to_HR2LQ.h"
#include "SEAPopularHRLQ.h"
#include "LpApproxSMFQ.h"
#include "Exact_Exponential_SMFQ.h"
#include "DirectApproachHR2LQ.h"
#include "Statistics.h"
#include "Popular.h"
#include "RHeuristicHRLQ.h"
#include "HHeuristicHRLQ.h"
#include "YokoiEnvyfreeHRLQ.h"
#include "MaximalEnvyfreeHRLQ.h"
#include "Utils.h"
#include "GraphReader.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <unistd.h>

std::stringstream stmp;

void compute_matching(bool A_proposing, const char* input_file, const char* output_file) {
    // setup input/output stream as std::cin/std::cout by default
    // if a file is specified use it to read/write

    auto cin_buf = std::cin.rdbuf(); // save pointer to std::cin buffer
    auto cout_buf = std::cout.rdbuf(); // save pointer to std::cout buffer

    std::ifstream filein(input_file);
    std::ofstream fileout(output_file);

    if (input_file) {
        std::cin.rdbuf(filein.rdbuf());
    }

    if (output_file) {
        std::cout.rdbuf(fileout.rdbuf());
    }

    std::shared_ptr<BipartiteGraph> G = GraphReader(std::cin).read_graph();
    if (G == NULL) {
        return;
    }
    MaxCardPopular alg(G, A_proposing);
    auto M = alg.compute_matching();

    // To get statistics of output matching
    // uncomment below lines

    //Statistics s;
    //s.get_statistics(G, M);
    //s.get_smfq_statistics(G, M);

    alg.checker(G, M, A_proposing,std::cout);

    // restore buffers
    std::cin.rdbuf(cin_buf);
    std::cout.rdbuf(cout_buf);
}

int main(int argc, char* argv[]) {
    int c = 0;
    bool A_proposing = true;
    const char* input_file = nullptr;
    const char* output_file = nullptr;

    opterr = 0;
    // choose the proposing partition using -A and -B
    // -o is the path where the matching computed should be stored
    while ((c = getopt(argc, argv, "ABo:")) != -1) {
        switch (c) {
            case 'A': A_proposing = true; break;
            case 'B': A_proposing = false; break; 
            case 'o': output_file = optarg; break;
            case '?':
                if (optopt == 'o') {
                    std::cerr << "Option -o requires an argument.\n";
                } else {
                    std::cerr << "Unknown option: " << static_cast<char>(optopt) << '\n';
                }
                break;
            default: break;
        }
    }

    // std::ofstream filelog(output_file);    
    for(int i=1; i<=1; i++)
    {
        char inp_file[] = "../testsuite/input/OneToOneX/TCY.txt";
        inp_file[27] = '0' + i;
        for(int j=0; j<1; j++)
        {
            inp_file[31] = '0' + j;  
            input_file = inp_file;

            // stmp << input_file << " : ";
            compute_matching(A_proposing, input_file, output_file);
        }
    }

    // for(int i=1; i<=4; i++)
    // {
    //     char inp_file[] = "../testsuite/input/ManyToOneX/TCY.txt";
    //     inp_file[28] = '0' + i;
    //     for(int j=0; j<10; j++)
    //     {
    //         inp_file[32] = '0' + j;  
    //         input_file = inp_file;
            
    //         // stmp << input_file << " : ";
    //         compute_matching(A_proposing, input_file, output_file);
    //     }
    // }

    // filelog << stmp.str();

    return 0;
}
