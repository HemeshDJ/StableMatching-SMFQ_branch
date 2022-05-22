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
#include <unistd.h>

template<typename T>
void compute_matching(bool A_proposing, bool signature, const char* input_file, const char* output_file, const char* log_file = nullptr) {
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
    T alg(G, A_proposing);
    auto M = alg.compute_matching();

    // To get statistics of output matching
    // uncomment below lines

    //Statistics s;
    //s.get_statistics(G, M);
    //s.get_smfq_statistics(G, M);
    
    if(signature) {
        print_signature(G, M, std::cout);
    }
    else 
        print_matching(G, M, std::cout);

    // restore buffers
    std::cin.rdbuf(cin_buf);
    std::cout.rdbuf(cout_buf);
}

int main(int argc, char* argv[]) {
    int c = 0;
    bool compute_sea_popular = false;
    bool compute_direct_sm2lq = false;
    bool compute_lp_smfq = false;
    bool compute_exact_exp_smfq = false;
    bool compute_stable = false;
    bool compute_popular = false;
    bool compute_max_card = false;
    bool compute_rhrlq = false;
    bool compute_hhrlq = false;
    bool compute_yhrlq = false;
    bool compute_ehrlq = false;
    bool convert_hr_to_hr2lq = false;
    bool signature = false;
    bool A_proposing = true;
    bool run_test_suite = false;
    const char* input_file = nullptr;
    const char* output_file = nullptr;

    opterr = 0;
    // choose the proposing partition using -A and -B
    // -s, -p, and -m flags compute the stable, max-card popular and pop among
    // max-card matchings respectively
    // -r and -h compute the resident and hopsital heuristic for an HRLQ instance
    // -i is the path to the input graph, -o is the path where the matching
    // computed should be stored
    // -t is to run the test suite and return the tests that failed
    while ((c = getopt(argc, argv, "ABczdklspmrhyegti:o:")) != -1) {
        switch (c) {
            case 'A': A_proposing = true; break;
            case 'B': A_proposing = false; break; 
            case 'c': convert_hr_to_hr2lq = true; break;
            case 'd': compute_direct_sm2lq = true; break;
            case 'k': compute_exact_exp_smfq = true; break;
            case 'l': compute_lp_smfq = true; break;
            case 'z': compute_sea_popular = true; break;
            case 's': compute_stable = true; break;
            case 'p': compute_popular = true; break;
            case 'm': compute_max_card = true; break;
            case 'r': compute_rhrlq = true; break;
            case 'h': compute_hhrlq = true; break;
            case 'y': compute_yhrlq = true; break;
            case 'e': compute_ehrlq = true; break;
            case 'g': signature = true; break;
            case 't': run_test_suite = true; break;
            case 'i': input_file = optarg; break;
            case 'o': output_file = optarg; break;
            case '?':
                if (optopt == 'i' && !run_test_suite) {
                    std::cerr << "Option -i requires an argument.\n";
                } else if (optopt == 'o') {
                    std::cerr << "Option -o requires an argument.\n";
                } else {
                    std::cerr << "Unknown option: " << static_cast<char>(optopt) << '\n';
                }
                break;
            default: break;
        }
    }
    char matching_dump[] = "../resources/test/dumpfile.txt";
    for(int i=1; i<=4; i++)
    {
        char inp_file[] = "../testcases/OneToOneX/TCY.txt";
        inp_file[18] = '0' + i;
        for(int j=0; j<10; j++)
        {
            inp_file[22] = '0' + j;  
            input_file = inp_file;

            compute_matching<MaxCardPopular>(A_proposing, signature, input_file, output_file);
        }
    }

    return 0;
}
