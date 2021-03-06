#include "Utils.h"
#include "Vertex.h"
#include "GraphReader.h"
#include <set>
#include <sstream>

const char* token_to_string(Token tok) {
   static struct {
       const char* desc;
       Token token;
   } array [] = {
           {"@", TOK_AT},
           {"@PartitionA", TOK_PARTITION_A},
           {"@PartitionB", TOK_PARTITION_B},
           {"@PreferenceListsA", TOK_PREF_LISTS_A},
           {"@PreferenceListsB", TOK_PREF_LISTS_B},
           {"@End", TOK_END},
           {"STRING", TOK_STRING},
           {":", TOK_COLON},
           {",", TOK_COMMA},
           {";", TOK_SEMICOLON},
           {"(", TOK_LEFT_BRACE},
           {")", TOK_RIGHT_BRACE},
           {"EOF", TOK_EOF},
           {"ERROR", TOK_ERROR},
           {nullptr}};

   for (int i = 0; array[i].desc != nullptr; ++i) {
       if (array[i].token == tok) {
           return array[i].desc;
       }
   }

   return "UNKNOWN";
}

bool is_number(const std::string& s) {
    for (int i = 0; i < s.length(); i++)
        if (s[i] < '0' || s[i] > '9')return false;
    return true;
}

int to_integer(const std::string& s) {
    return (int) std::strtol(s.c_str(), nullptr, 10);
}

PartnerList get_partners(std::shared_ptr<MatchingAlgorithm::MatchedPairListType> M, VertexPtr v) {
    auto v_it = M->find(v);
    return v_it == M->end() ? PartnerList() : v_it->second;
}

PartnerList::SizeType number_of_partners(std::shared_ptr<MatchingAlgorithm::MatchedPairListType> M, VertexPtr v) {
    return get_partners(std::move(M), std::move(v)).size();
}

bool has_partner(std::shared_ptr<MatchingAlgorithm::MatchedPairListType> M, VertexPtr v) {
    return number_of_partners(std::move(M), std::move(v)) > 0;
}

VertexPtr get_partner(std::shared_ptr<MatchingAlgorithm::MatchedPairListType> M, VertexPtr v) {
    return get_partner(get_partners(std::move(M), std::move(v)));
}

VertexPtr get_partner(const PartnerList& partner_list) {
    if (partner_list.empty()) {
        return nullptr;
    } else {
        return partner_list.get_least_preferred().vertex;
    }
}

int matching_size(std::shared_ptr<MatchingAlgorithm::MatchedPairListType> M) {
    int size = 0;

    for (auto& it : *M) {
        size += it.second.size();
    }

    return size / 2;
}

VertexPtr get_vertex_by_id(std::shared_ptr<BipartiteGraph> G, const IdType& id) {
    auto A_it = G->get_A_partition().find(id);

    return A_it == G->get_A_partition().end() ? G->get_B_partition().find(id)->second : A_it->second;
}

std::shared_ptr<BipartiteGraph> read_graph(const std::string& filepath) {
    std::ifstream fin(filepath);
    return GraphReader(fin).read_graph();
}

// a new id is of the form id^k
IdType get_vertex_id(const IdType& id, int k) {
    std::stringstream stmp;
    stmp << id << '^' << k;
    return stmp.str();
}

int get_vertex_level(const IdType& id) {
    return to_integer(id.substr(id.find('^')+1));
}

// a dummy id is of the form d^k_id
IdType get_dummy_id(const IdType& id, int k) {
    std::stringstream stmp;
    stmp << "d^" << k << '_' << id;
    return stmp.str();
}

// return the dummy level from the given id
int get_dummy_level(const IdType& id) {
    auto caret_pos = id.find('^');
    auto underscore_pos = id.find('_');
    return to_integer(id.substr(caret_pos+1, underscore_pos-caret_pos-1));
}

RankType compute_rank(VertexPtr u, const PreferenceList& pref_list) {
    auto index = pref_list.find_index(std::move(u));
    return (RankType) ((index + 1) /*+ (level * pref_list.size())*/ );
}

void print_matching(std::shared_ptr<BipartiteGraph> G,
                    std::shared_ptr<MatchingAlgorithm::MatchedPairListType> M, std::ostream& out)
{
    std::stringstream stmp;

    for (const auto& it : G->get_A_partition()) {
        auto u = it.second;
        auto M_u = M->find(u);

        if (M_u != M->end()) {
            auto& partners = M_u->second;

            for (const auto& i : partners) {
                auto v = i.vertex;

                stmp << u->get_id() << ','
                     << v->get_id() << ','
                     << i.rank << '\n';
            }
        }
    }

    out << stmp.str();
}

void print_signature(std::shared_ptr<BipartiteGraph> G,
                    std::shared_ptr<MatchingAlgorithm::MatchedPairListType> M, std::ostream& out)
{
    std::stringstream stmp;

    int unmatched = (G->get_A_partition()).size();
    int num_ranks = 1;

    for (const auto& it : G->get_A_partition()) {
        auto u = it.second;
        if(num_ranks < u->get_preference_list().size()+1) 
            num_ranks = u->get_preference_list().size()+1;
    }

    int ranks[num_ranks];
    for(int i = 0; i < num_ranks; i++) {
        ranks[i] = 0;
    }

    for (const auto& it : G->get_A_partition()) {
        auto u = it.second;
        auto M_u = M->find(u);

        if (M_u != M->end()) {
            auto& partners = M_u->second;

            for (const auto& i : partners) {
                ranks[i.rank]++;
                unmatched--;
            }
        }
    }
    stmp << "Signature" << "\n";
    stmp << "|A| = " << (G->get_A_partition()).size() << "\n";
    stmp << "|B| = " << (G->get_B_partition()).size() << "\n";
    stmp << "Max Pref Size of A = " << num_ranks-1 << "\n";
    for(int i = 1; i < num_ranks; i++) {
        stmp << i << " " << ranks[i] << "\n";
    }
    
    // Instead of last rank
    stmp << "Unassigned " << unmatched << "\n";

    out << stmp.str();
}   