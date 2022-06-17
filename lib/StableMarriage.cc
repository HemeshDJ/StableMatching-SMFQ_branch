#include "StableMarriage.h"
#include "Vertex.h"
#include "Partner.h"
#include "Utils.h"
#include <sstream>

void StableMarriage::checker(std::shared_ptr<BipartiteGraph> G,std::shared_ptr<MatchingAlgorithm::MatchedPairListType> M, bool A_proposing, std::ostream &out)
{
    std::stringstream stmp;

    auto proposing_partition = A_proposing ? G->get_A_partition() 
                                                : G->get_B_partition();

    // assigning edge weights to the edges
    std::map<VertexPtr, std::map<VertexPtr, int>> edge_weights;

    for( auto & it : proposing_partition ) {
        auto u = it.second;
        auto u_pref_list = u->get_preference_list();

        for( auto& it2 : u_pref_list ) {
            int weight = 0;
            auto v = it2.vertex;
            
            auto v_pref_list = v->get_preference_list();
            auto M_v = M->find(v);
            if(M_v != M->end()) {
                auto& partners = M_v->second;
                if(partners.size() == 0) {
                    weight += 1;
                }
                else {
                    auto v_worst_partner = partners.get_least_preferred();
                    auto u_rank = compute_rank(u, v_pref_list);
                    
                    if(v_worst_partner.rank > u_rank) {
                        weight += 1;
                    }
                    else if(u_rank > v_worst_partner.rank) {
                        weight += -1;
                    }
                    else {
                        weight += 0;
                    }
                }
            }
            else {
                weight += 1;
            }

            auto M_u = M->find(u);
            if(M_u != M->end()) {
                auto& partners = M_u->second;
                if(partners.size() == 0) {
                    weight += 1;
                }
                else {
                    auto u_worst_partner = partners.get_least_preferred();
                    auto v_rank = compute_rank(v, u_pref_list);
                    
                    if(u_worst_partner.rank > v_rank) {
                        weight += 1;
                    }
                    else if(v_rank > u_worst_partner.rank) {
                        weight += -1;
                    }
                    else {
                        weight += 0;
                    }
                }
            }
            else {
                weight += 1;
            }
            
            edge_weights[u][v] += weight;
            edge_weights[v][u] += weight;
        }
    }

    // checking if each edge is covered
    int flag = 1;
    for( auto & it : proposing_partition ) {
        auto u = it.second;

        for( auto& it2 : u->get_preference_list() ) {
            auto v = it2.vertex;

            if(edge_weights[u][v] == 2 && v->get_upper_quota()) {           //if u prefers v over M(u) and vice-versa
                
                if(flag)
                    stmp << "Matching is not stable" << std::endl;
                
                flag = 0;
                stmp << "Edge " << u->get_id() << " -- " << v->get_id() << " is a blocking pair!" << std::endl;
            }
        }
    }


    // checking if total sum of dual_value is 0
    if(flag) {
        stmp << "Matching is Stable" << std::endl;
    }

    out << stmp.str();

}