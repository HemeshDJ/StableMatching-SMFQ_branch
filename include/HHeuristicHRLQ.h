#ifndef HHEURISTIC_HRLQ_H
#define HHEURISTIC_HRLQ_H

#include "MatchingAlgorithm.h"

class HHeuristicHRLQ : public MatchingAlgorithm {
private:
    // graph for phase 2
    std::shared_ptr<BipartiteGraph> G2_;

    // matching to hold temporary results
    std::shared_ptr<MatchedPairListType> M_tmp_;

    bool compute_phase2_matching(std::shared_ptr<MatchedPairListType> M,
                                 std::map<IdType, PartnerList::SizeType >& def,
                                 int s);
    // augment graph phase 2
    std::shared_ptr<BipartiteGraph> augment_phase2(int s);

public:
    explicit HHeuristicHRLQ(std::shared_ptr<BipartiteGraph> G, bool A_proposing=true);
    ~HHeuristicHRLQ() override = default;

    std::shared_ptr<MatchedPairListType> compute_matching() override;

    void check_popularity(std::shared_ptr<BipartiteGraph> G, 
        std::shared_ptr<MatchingAlgorithm::MatchedPairListType> M, bool A_proposing, std::ostream& out) 
    {}
};

#endif
