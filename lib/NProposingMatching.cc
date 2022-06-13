#include "NProposingMatching.h"
#include "Vertex.h"
#include "Partner.h"
#include "Utils.h"
#include <set>
#include <sstream>

NProposingMatching::NProposingMatching(std::shared_ptr<BipartiteGraph> G,
                                       bool A_proposing, int max_level)
    : MatchingAlgorithm(std::move(G), A_proposing), max_level(max_level)
{}

VertexPtr NProposingMatching::remove_and_update(FreeListType& free_list,
                                                std::map<VertexPtr,
                                                VertexBookkeeping>& bookkeep_data) {
    auto ret = free_list.top();
    free_list.pop(); // remove u from free_list
    bookkeep_data[ret].in_free_list = false;
    return ret;
}

void NProposingMatching::add_to_list(FreeListType& free_list, VertexPtr u) {
    free_list.push(std::move(u));
}

void NProposingMatching::add_to_list_and_update(NProposingMatching::FreeListType &free_list,
                                                VertexBookkeeping &u_data, VertexPtr u) {
    if (not u_data.in_free_list) {
        u_data.begin += 1;
        u_data.in_free_list = true;
        add_to_list(free_list, std::move(u));
    }
}

void NProposingMatching::add_matched_partners(std::shared_ptr<MatchedPairListType> M,
                                              VertexPtr u, VertexPtr v,
                                              const VertexBookkeeping& u_data,
                                              const PreferenceList& v_pref_list) {
    // v is the first vertex on u's current preference list
    // invariant for non proposing vertices: rank = index + 1
    // v is always at level 0
    add_partner(M, u, v, (RankType) u_data.begin + 1, 0);
    add_partner(M, v, u, compute_rank(u, v_pref_list), u_data.level);
}

void NProposingMatching::checker(std::shared_ptr<BipartiteGraph> G,std::shared_ptr<MatchingAlgorithm::MatchedPairListType> M, bool A_proposing, std::ostream& out)
{
    std::stringstream stmp;
    const auto& proposing_partition = A_proposing ? G->get_A_partition()
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

    for (const auto& it : proposing_partition) {
        auto u = it.second;
        auto M_u = M->find(u);

        if (M_u != M->end()) {
            auto& partners = M_u->second;

            for (const auto& i : partners) {
                auto v = i.vertex;
                edge_weights[v][v] = -1;
                edge_weights[u][u] = -1;
            }
        }
    }

    // assigning popularity to each vertex
    std::map<VertexPtr, int> dual_value;    

    for (const auto& it : proposing_partition) {
        auto u = it.second;
        auto M_u = M->find(u);

        if (M_u != M->end()) {
            auto& partners = M_u->second;

            for (const auto& i : partners) {
                auto v = i.vertex;
                dual_value[u] = (levels[u] ? -1 : 1);
                dual_value[v] += (levels[u] ? 1 : -1);
            }

            if(partners.size() == 0) {
                dual_value[u] = 0;
            }
        }
        else {
            dual_value[u] = 0;
        }
    }                                      

    // checking if each edge is covered
    int flag = 1;
    for( auto & it : proposing_partition ) {
        auto u = it.second;

        for( auto& it2 : u->get_preference_list() ) {
            auto v = it2.vertex;

            if(dual_value[u] + dual_value[v] < edge_weights[u][v]) {
                flag = 0;
                stmp << u->get_id() << "," << v->get_id() << std::endl;
            }
        }
    }


    // checking if total sum of dual_value is 0
    int dual_sum = 0;
    for(auto it : dual_value) {
        dual_sum += it.second;
    }
    
    if(!flag)
        stmp << "Edges not covered. " << std::endl;
    if(dual_sum)
        stmp << "Popularity sum is non-zero, equal to " << dual_sum << ". " << std::endl;
    if(flag && !dual_sum)
        stmp << "Passed. Certificate issued!" << std::endl;

    out << stmp.str();
}

void NProposingMatching::checker(std::shared_ptr<BipartiteGraph> G,std::shared_ptr<MatchingAlgorithm::MatchedPairListType> M, bool A_proposing, std::ostream& out)
{
    std::stringstream stmp;
    const auto& proposing_partition = A_proposing ? G->get_A_partition()
                                                       : G->get_B_partition();

    BipartiteGraph::ContainerType Am, Bm;

    for(auto &it : G->get_A_partition()) {
        auto u = it.second;
        Am.emplace(u->get_id(), std::make_shared<Vertex>(u->get_id(), 0, 1));
    }

    std::map<VertexPtr,int> availableB;
    for(auto &it : G->get_B_partition()) {
        auto u = it.second;
        availableB[u] = 0;
        for(int i = 0; i < u->get_upper_quota(); i++) {
            std::string u_id = u->get_id() + "///" + std::to_string(i);
            Bm.emplace(u_id, std::make_shared<Vertex>(u_id, 0, 1));
        }
    }

    for(auto &it : G->get_A_partition()) {
        auto u = it.second;
        auto u_pref_list = u->get_preference_list();
        auto M_u = M->find(u);
        PreferenceList& pref_list = Am[u]->get_preference_list();

        if(M_u != M->end()) {
            auto& partners = M_u->second;
            for(auto &it2 : u_pref_list) {
                auto v = it2.vertex;
                if(partners.find(v) != partners.end())
                {
                    std::string v_id = v->get_id() + "///" + std::to_string(availableB[v]);
                    pref_list.emplace_back(Bm[v_id]);
                    (Bm[v_id]->get_preference_list()).emplace_back(u.rank, u);
                    availableB[v]++;
                }
                else 
                {
                    for(int i = 0; i < v->get_upper_quota(); i++) {
                        std::string v_id = v->get_id() + "///" + std::to_string(i);
                        pref_list.emplace_back(Bm[v_id]);
                        (Bm[v_id]->get_preference_list()).emplace_back(u.rank, u);
                    }
                }   
            }
        }
        else {
            for(auto &it2 : u_pref_list) {
                auto v = it2.vertex;
                for(int i = 0; i < v->get_upper_quota(); i++) {
                    std::string v_id = v->get_id() + "///" + std::to_string(i);
                    pref_list.emplace_back(Bm[v_id]);
                }
            }
        }
    }

    // for(auto &it : G->get_B_partition()) {
    //     auto u = it.second;
    //     auto& u2 = u->get_id();
    //     while(u2.back() != '/')
    //         u2.pop_back();
    //     u2.pop_back();u2.pop_back();u2.pop_back();
        
    //     auto u_pref_list = u2->get_preference_list();
    //     auto M_u = M->find(u);

    //     if(M_u != M->end()) {
    //         auto& partners = M_u->second;
    //         for(auto &it2 : u_pref_list) {
    //             auto v = it2.vertex;
    //             while(v->get_id().back() != '/')
    //                 v->get_id()->pop_back();
    //             v->get_id()->pop_back();
    //             v->get_id()->pop_back();
    //             v->get_id()->pop_back();
    //             if(partners.find(v) != partners.end())
    //             {
    //                 PreferenceList& vm_pref_list = Am[v]->get_preference_list();
    //                 if(vm_pref_list.find(u) != vm_pref_list.end()) {
    //                     pref_list.emplace_back(v);
    //                 }
    //             }       
    //             else 
    //             {
    //                 pref_list.emplace_back(v);
    //             } 
    //         }
    //     }
    //     else {
    //         for(auto &it2 : u_pref_list) {
    //             auto v = it2.vertex;
    //             for(int i = 0; i < v->get_upper_quota(); i++) {
    //                 std::string v_id = v->get_id() + "///" + std::to_string(i);
    //                 pref_list.emplace_back(Am[v_id]);
    //             }
    //         }
    //     }
    // }
    
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

    for (const auto& it : proposing_partition) {
        auto u = it.second;
        auto M_u = M->find(u);

        if (M_u != M->end()) {
            auto& partners = M_u->second;

            for (const auto& i : partners) {
                auto v = i.vertex;
                edge_weights[v][v] = -1;
                edge_weights[u][u] = -1;
            }
        }
    }

    // assigning popularity to each vertex
    std::map<VertexPtr, int> dual_value;    

    for (const auto& it : proposing_partition) {
        auto u = it.second;
        auto M_u = M->find(u);

        if (M_u != M->end()) {
            auto& partners = M_u->second;

            for (const auto& i : partners) {
                auto v = i.vertex;
                dual_value[u] = (levels[u] ? -1 : 1);
                dual_value[v] += (levels[u] ? 1 : -1);
            }

            if(partners.size() == 0) {
                dual_value[u] = 0;
            }
        }
        else {
            dual_value[u] = 0;
        }
    }                                      

    // checking if each edge is covered
    int flag = 1;
    for( auto & it : proposing_partition ) {
        auto u = it.second;

        for( auto& it2 : u->get_preference_list() ) {
            auto v = it2.vertex;

            if(dual_value[u] + dual_value[v] < edge_weights[u][v]) {
                flag = 0;
                stmp << u->get_id() << "," << v->get_id() << std::endl;
            }
        }
    }


    // checking if total sum of dual_value is 0
    int dual_sum = 0;
    for(auto it : dual_value) {
        dual_sum += it.second;
    }
    
    if(!flag)
        stmp << "Edges not covered. " << std::endl;
    if(dual_sum)
        stmp << "Popularity sum is non-zero, equal to " << dual_sum << ". " << std::endl;
    if(flag && !dual_sum)
        stmp << "Passed. Certificate issued!" << std::endl;

    out << stmp.str();
}

std::shared_ptr<MatchingAlgorithm::MatchedPairListType> NProposingMatching::compute_matching() {
    FreeListType free_list;
    std::map<VertexPtr, VertexBookkeeping> bookkeep_data;
    std::shared_ptr<BipartiteGraph> G = get_graph();
    auto M = std::make_shared<MatchingAlgorithm::MatchedPairListType>();

    // choose the partitions from which the vertices will propose
    const auto& proposing_partition = is_A_proposing() ? G->get_A_partition()
                                                       : G->get_B_partition();

    // set the level of every vertex in the proposing partition to 0
    // mark all proposing vertices free (by pushing into the free_list)
    // and vertices from the opposite partition implicitly free
    for (auto& it : proposing_partition) {
        auto v = it.second;
        free_list.push(v);
        bookkeep_data[v] = VertexBookkeeping(0, v->get_preference_list().size());
    }

    // there is at least one vertex in the free list
    while (not free_list.empty()) {
        // arbitrary vertex in free list
        auto u = remove_and_update(free_list, bookkeep_data);
        auto& u_pref_list = u->get_preference_list();

        // if u^l can have a partner and hasn't exhausted its preference list
        if (u->get_upper_quota() > 0 and not bookkeep_data[u].is_exhausted()) {
            // highest ranked vertex to whom u has not yet proposed
            auto v = u_pref_list.at(bookkeep_data[u].begin).vertex;
            auto v_pref_list = v->get_preference_list();

            
            // v cannot be matched to anyone
            if (v->get_upper_quota() == 0 or v_pref_list.size() == 0) {
                // do nothing, also inconsistent graph maybe
            } else if (number_of_partners(M, v) == v->get_upper_quota()) {
                // |M[v]| = upper_quota(v)
                auto v_worst_partner = M->at(v).get_least_preferred();
                auto possible_partner = Partner(u, compute_rank(u, v_pref_list), bookkeep_data[u].level);
                
                if (v_worst_partner < possible_partner) {
                    // remove M[v_worst_partner] from M[v], and M[v] from M[v_worst_partner]
                    M->at(v).remove_least_preferred();
                    M->at(v_worst_partner.vertex).remove(v);

                    // add u and v to the matching
                    add_matched_partners(M, u, v, bookkeep_data[u], v_pref_list);

                    // add v_worst_partner to free_list
                    add_to_list_and_update(free_list,
                                           bookkeep_data[v_worst_partner.vertex],
                                           v_worst_partner.vertex);
                }
            } else {
                // add u and v to the matching
                add_matched_partners(M, u, v, bookkeep_data[u], v_pref_list);
            }

            // add u to the free_list if it has residual capacity
            if (u->get_upper_quota() > number_of_partners(M, u)) {
                add_to_list_and_update(free_list, bookkeep_data[u], u);
            }
        } else if (bookkeep_data[u].level < max_level) {
            bookkeep_data[u].level += 1;
            bookkeep_data[u].begin = 0; // reset proposal index
            bookkeep_data[u].in_free_list = true;
            add_to_list(free_list, u);
        }
    }
    
    /** New Code **/

    for (const auto& it : proposing_partition) {
        auto u = it.second;
        auto M_u = M->find(u);

        if (M_u != M->end()) {
            auto& partners = M_u->second;

            for (const auto& i : partners) {
                auto v = i.vertex;
                levels[u] = bookkeep_data[u].level;
            }

            if (partners.size() == 0) {
                levels[u] = bookkeep_data[u].level;
            }
        }
        else {
            levels[u] = bookkeep_data[u].level;
        }
    }

    /** End New Code **/

    return M;
}