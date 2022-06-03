#ifndef POPULAR_H
#define POPULAR_H

#include "MatchingAlgorithm.h"
#include "NProposingMatching.h"

/// max-card popular matching in an HR instance
class MaxCardPopular : public MatchingAlgorithm {
public:
    NProposingMatching npm = NProposingMatching(get_graph(), is_A_proposing(), 1);
    
    explicit MaxCardPopular(std::shared_ptr<BipartiteGraph> G,
                   bool A_proposing=true)
        : MatchingAlgorithm(G, A_proposing)
    {}

    ~MaxCardPopular() override = default;

    std::shared_ptr<MatchedPairListType> compute_matching() override {
        // NProposingMatching npm (get_graph(), is_A_proposing(), 1);
        return npm.compute_matching();
    }

    std::pair<int,int> checker(std::shared_ptr<BipartiteGraph> G, 
        std::shared_ptr<MatchingAlgorithm::MatchedPairListType> M, bool A_proposing) {
        // NProposingMatching npm (get_graph(), is_A_proposing(), 1);
        return npm.checker(G, M, A_proposing);
    }

};

/// a popular matching among the set of max-card matchings in an instance
class PopularAmongMaxCard  : public MatchingAlgorithm {
public:
    explicit PopularAmongMaxCard(std::shared_ptr<BipartiteGraph> G,
                   bool A_proposing=true)
        : MatchingAlgorithm(G, A_proposing)
    {}

    ~PopularAmongMaxCard() override = default;

    std::shared_ptr<MatchedPairListType> compute_matching() override {
        NProposingMatching npm (get_graph(), is_A_proposing(), 10);//(int) get_graph()->get_A_partition().size() - 1);
        return npm.compute_matching();
    }

    std::pair<int,int> checker(std::shared_ptr<BipartiteGraph> G, 
        std::shared_ptr<MatchingAlgorithm::MatchedPairListType> M, bool A_proposing) 
    { return {-1,-1}; }
};

#endif
