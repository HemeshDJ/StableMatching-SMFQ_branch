#ifndef RHEURISTIC_HRLQ_H
#define RHEURISTIC_HRLQ_H

#include "MatchingAlgorithm.h"

class RHeuristicHRLQ : public MatchingAlgorithm {
private:
    // graphs for phase 1 and 2
    std::shared_ptr<BipartiteGraph> G1_;
    std::shared_ptr<BipartiteGraph> G2_;

    // augment graph for phase 1 and phase 2
    std::shared_ptr<BipartiteGraph> augment_phase1();
    std::shared_ptr<BipartiteGraph> augment_phase2(std::shared_ptr<MatchedPairListType> M);

public:
    // A_proposing should always be true, this is a resident proposing
    // heuristic
    explicit RHeuristicHRLQ(std::shared_ptr<BipartiteGraph> G, bool A_proposing=true);
    ~RHeuristicHRLQ() override = default;

    std::shared_ptr<MatchedPairListType> compute_matching() override;

    std::pair<int,int> check_popularity(std::shared_ptr<BipartiteGraph> G, 
        std::shared_ptr<MatchingAlgorithm::MatchedPairListType> M, bool A_proposing, std::ostream& out) 
    { return {-1,-1}; }
};

#endif
