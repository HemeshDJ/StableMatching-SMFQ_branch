#include "BipartiteGraph.h"
#include "StableMarriage.h"
// #include "Convert_HR_to_HR2LQ.h"
// #include "SEAPopularHRLQ.h"
// #include "LpApproxSMFQ.h"
// #include "Exact_Exponential_SMFQ.h"
// #include "DirectApproachHR2LQ.h"
#include "Statistics.h"
#include "Popular.h"
// #include "RHeuristicHRLQ.h"
// #include "HHeuristicHRLQ.h"
// #include "YokoiEnvyfreeHRLQ.h"
// #include "MaximalEnvyfreeHRLQ.h"
#include "Utils.h"
#include "GraphReader.h"
#include <fstream>
#include <iostream>
#include <unistd.h>

template<typename T>
void compute_matching(bool A_proposing, bool test, const char* input_file, const char* output_file, const char* sig_file = nullptr) {
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
    
    if(test) {
        alg.checker(G, M, A_proposing, std::cerr);
    }

    print_matching(G, M, std::cout);

    // restore buffers
    std::cin.rdbuf(cin_buf);
    std::cout.rdbuf(cout_buf);
    filein.close();
    fileout.close();

    if(sig_file) {
        fileout.open(sig_file);
        std::cout.rdbuf(fileout.rdbuf());
        print_signature(G, M, std::cout);
        std::cout.rdbuf(cout_buf);
        fileout.close();
    }
}

int main(int argc, char* argv[]) {
    int c = 0;
    // bool compute_sea_popular = false;
    // bool compute_direct_sm2lq = false;
    // bool compute_lp_smfq = false;
    // bool compute_exact_exp_smfq = false;
    bool compute_stable = false;
    bool compute_popular = false;
    bool compute_max_card = false;
    // bool compute_rhrlq = false;
    // bool compute_hhrlq = false;
    // bool compute_yhrlq = false;
    // bool compute_ehrlq = false;
    // bool convert_hr_to_hr2lq = false;
    bool test = false;
    bool A_proposing = true;
    const char* input_file = nullptr;
    const char* output_file = nullptr;
    const char* sig_file = nullptr;

    opterr = 0;
    // choose the proposing partition using -A and -B
    // -s, -p, and -m flags compute the stable, max-card popular and pop among
    // max-card matchings respectively
    // -r and -h compute the resident and hopsital heuristic for an HRLQ instance
    // -g prints the signature of the matching
    // -t computes in test mode
    // -i is the path to the input graph, -o is the path where the matching
    // computed should be stored
    while ((c = getopt(argc, argv, "ABspmtg:i:o:")) != -1) {
        switch (c) {
            case 'A': A_proposing = true; break;
            case 'B': A_proposing = false; break; 
            case 's': compute_stable = true; break;
            case 'p': compute_popular = true; break;
            case 't': test = true; break;
            // case 'c': convert_hr_to_hr2lq = true; break;
            // case 'd': compute_direct_sm2lq = true; break;
            // case 'k': compute_exact_exp_smfq = true; break;
            // case 'l': compute_lp_smfq = true; break;
            // case 'z': compute_sea_popular = true; break;
            case 'm': compute_max_card = true; break;
            // case 'r': compute_rhrlq = true; break;
            // case 'h': compute_hhrlq = true; break;
            // case 'y': compute_yhrlq = true; break;
            // case 'e': compute_ehrlq = true; break;
            case 'g': sig_file = optarg; break;
            case 'i': input_file = optarg; break;
            case 'o': output_file = optarg; break;
            case '?':
                if (optopt == 'i') {
                    std::cerr << "Option -i requires an argument.\n";
                } else if (optopt == 'o') {
                    std::cerr << "Option -o requires an argument.\n";
                } else if (optopt == 'g') {
                    std::cerr << "Option -g requires an argument.\n";
                } else {
                    std::cerr << "Unknown option: " << static_cast<char>(optopt) << '\n';
                }
                break;
            default: break;
        }
    }

    if(!input_file) {
        std::cerr << "Add -i <input_file_name> to specify input file.\n";
        return 0;
    }
    if(!output_file && !sig_file) {     //neither signature nor output file is specified
        std::cerr << "Add -o or -g (with arguments).\n";
        return 0;
    }
    
    if (compute_stable) {
        compute_matching<StableMarriage>(A_proposing, test, input_file, output_file, sig_file);
    // }else if (convert_hr_to_hr2lq) {
    //     compute_matching<Convert_HR_to_HR2LQ>(A_proposing, test, input_file, output_file);
    // }else if (compute_direct_sm2lq) {
    //     compute_matching<DirectApproachHR2LQ>(A_proposing, test, input_file, output_file);
    // }else if (compute_exact_exp_smfq) {
    //     compute_matching<Exact_Exponential_SMFQ>(A_proposing, test, input_file, output_file);
    // }else if (compute_lp_smfq) {
    //     compute_matching<LpApproxSMFQ>(A_proposing, test, input_file, output_file);
    // }else if (compute_sea_popular) {
    //     compute_matching<SEAPopularHRLQ>(A_proposing, test, input_file, output_file);
    }else if (compute_popular) {
        compute_matching<MaxCardPopular>(A_proposing, test, input_file, output_file, sig_file);
    } else if (compute_max_card) {
        compute_matching<PopularAmongMaxCard>(A_proposing, test, input_file, output_file);
    // } else if (compute_rhrlq) {
    //     compute_matching<RHeuristicHRLQ>(A_proposing, test, input_file, output_file);
    // } else if (compute_hhrlq) {
    //     compute_matching<HHeuristicHRLQ>(A_proposing, test, input_file, output_file);
    // } else if (compute_yhrlq) {
    //     compute_matching<YokoiEnvyfreeHRLQ>(A_proposing, test, input_file, output_file);
    // } else if (compute_ehrlq) {
    //     compute_matching<MaximalEnvyfreeHRLQ>(A_proposing, test, input_file, output_file);
    }
    else {
        std::cerr << "Add -s to compute stable matching, -p for maximum cardinality popular matching. \n";
    }

    return 0;
}
