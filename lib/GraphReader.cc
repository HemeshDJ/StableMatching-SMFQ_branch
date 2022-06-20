#include "GraphReader.h"
#include "Vertex.h"
#include "PreferenceList.h"
#include "Utils.h"
#include <cstdlib>
#include <iostream>
#include <stdexcept>

/// Lexer class defined here

Lexer::Lexer(std::istream& in)
    : ch_(' '), lineno_(1), in_(in)
{
    if (not in_) {
        throw ReaderException("error reading input stream.");
    }
}

void Lexer::read_character() {
    if (ch_ == '\n') {
        ++lineno_;
    }
    ch_ = in_.get();
}

int Lexer::line_number() {
    return lineno_;
}

Token Lexer::next_token() {
    // skip whitespace
    while (isspace(ch_)) {
        read_character();
    }

    // skip comments
    if (ch_ == '#') {
        while (ch_ != '\n' and ch_ != EOF) {
            read_character();
        }
    }

    // skip newline and spaces
    if (isspace(ch_)) { read_character(); return next_token(); }
    if (ch_ == EOF)   { return TOK_EOF; }
    if (ch_ == ':')   { read_character(); return TOK_COLON; }
    if (ch_ == '@')   { read_character(); return TOK_AT; }
    if (ch_ == ',')   { read_character(); return TOK_COMMA; }
    if (ch_ == ';')   { read_character(); return TOK_SEMICOLON; }
    if (ch_ == '(')   { read_character(); return TOK_LEFT_BRACE; }
    if (ch_ == ')')   { read_character(); return TOK_RIGHT_BRACE; }

    // a directive or a string
    if (isalnum(ch_)) {
        lexeme_.clear();

        while (isalnum(ch_) or ch_ == '+') {
            lexeme_.push_back(static_cast<char>(ch_));
            read_character();
        }
        if (lexeme_ == "End") return TOK_END;
        if (lexeme_ == "PartitionA") return TOK_PARTITION_A;
        if (lexeme_ == "PartitionB") return TOK_PARTITION_B;
        if (lexeme_ == "PreferenceListsA") return TOK_PREF_LISTS_A;
        if (lexeme_ == "PreferenceListsB") return TOK_PREF_LISTS_B;
        return TOK_STRING;
    }

    // flag error, return the erraneous character
    lexeme_.clear();
    lexeme_.push_back(static_cast<char>(ch_));
    read_character();
    return TOK_ERROR;
}

std::string const& Lexer::get_lexeme() const {
    return lexeme_;
}

/// GraphReader class defined here

GraphReader::GraphReader(std::istream& in) {
    lexer_ = std::make_unique<Lexer>(in);
    consume(); // read first token
}

void GraphReader::consume() {
    curtok_ = lexer_->next_token();
}

void GraphReader::match(Token expected) {
    if (curtok_ != expected) {
        error_occurred = true;
        std::cout << "Line " << error_message("invalid data in file", curtok_, {expected}) << "\n";
        exit(0);
    //    consume();
    //    throw ReaderException(error_message("invalid data in file", curtok_, {expected}));
    } else {
        consume(); // skip the token
    }
}

const char* GraphReader::error_message(const char* prefix, Token got, const std::vector<Token>& expected) {
    static char buf[1024];
    int ret = snprintf(buf, sizeof buf, "%d: %s", lexer_->line_number(), prefix);
    ret += snprintf(buf + ret, sizeof buf, ", got: '%s', expected one of: {", token_to_string(got));

    for (std::vector<Token>::size_type i = 0; i < expected.size(); ++i) {
        // token description
        ret += snprintf(buf + ret, sizeof buf, "'%s'", token_to_string(expected[i]));
        // pretty printing
        ret += snprintf(buf + ret, sizeof buf, "%s", (i == expected.size() - 1) ? "}" : ", ");
    }

    return buf;
}

/// partition are of the format
/// @Partition (A|B)
/// a, b, c ;
/// @End
void GraphReader::read_partition(BipartiteGraph::ContainerType& vmap) {
    // read the vertices in the partion
    while (curtok_ != TOK_SEMICOLON) {
        std::string v = lexer_->get_lexeme();

        // if there is already a vertex with id v
        if (vmap.find(v) != vmap.end()) {
            error_occurred = true;
            std::cout <<"Line " <<lexer_->line_number() <<": Duplicate vertex : "<<v<<"\n";
        }

        int lower_quota = 0, upper_quota = 1;
        match(TOK_STRING);

        // does this vertex specify quota(s)
        // (upper) or (lower, upper)
        if (curtok_ == TOK_LEFT_BRACE) {
            // eat '('
            match(TOK_LEFT_BRACE);

            // if quota is not given as positive integer or 0
            if (!is_number(lexer_->get_lexeme())) {
                error_occurred = true;
                std::cout << "Line " << lexer_->line_number() << ": Expected number for quota for vertex : " << v << "\n";
            }

            // read the upper quota
            upper_quota = to_integer(lexer_->get_lexeme());
            match(TOK_STRING);

            // check if this vertex has a lower quota as well
            if (curtok_ == TOK_COMMA) {
                match(TOK_COMMA);

                // the quota read first was the lower quota
                lower_quota = upper_quota;

                // if quota is not given as positive integer or 0
                if (!is_number(lexer_->get_lexeme())) {
                    std::cout << "Line " << lexer_->line_number() << ": Expected number for quota for vertex : " << v << "\n";
                }

                upper_quota = to_integer(lexer_->get_lexeme());
                match(TOK_STRING);
            }

            // eat ')'
            match(TOK_RIGHT_BRACE);
        }

        // if lower quota is greater than upper quota
        if (lower_quota > upper_quota) {
            error_occurred = true;
            std::cout << "Line " << lexer_->line_number() << ": Lower quota cannot greater than Upper quota for vertex " << v << "\n";
        }

        // add this vertex with the required quotas
        vmap.emplace(v, std::make_shared<Vertex>(v, lower_quota, upper_quota));

        // if there are more vertices, they must
        // be delimited using commas
        if(curtok_ == TOK_AT) {
            match(TOK_SEMICOLON);
        }
        else if (curtok_ != TOK_SEMICOLON) {
            match(TOK_COMMA);
        }
    }

    // list should be delimited by a semicolon
    match(TOK_SEMICOLON);

    // end of this directive
    match(TOK_AT);
    match(TOK_END);
}

/// @PreferenceLists (A|B)
/// preference lists
/// @End
/// preference lists for a vertex are given in this format
/// v: a, b, c ;
/// we read the preference lists for vertices in partition A
/// and lookup the other side in partition B
void GraphReader::read_preference_lists(BipartiteGraph::ContainerType& A, BipartiteGraph::ContainerType& B) {
    // read the lists
    while (curtok_ != TOK_AT) {
        // read the vertex for which the preference list is given
        std::string a = lexer_->get_lexeme();

        // if there is no vertex with that id 
        if (A.find(a) == A.end()) {
            error_occurred = true;
            std::cout << "Line " << lexer_->line_number() << ": Vertex not found: " << a << "\n";

            match(TOK_STRING);
            match(TOK_COLON);

            while (curtok_ != TOK_SEMICOLON && curtok_ != TOK_EOF) {
                consume();
            }
            if (curtok_ == TOK_EOF)return;
            match(TOK_SEMICOLON);
            continue;
        }

        match(TOK_STRING);
        match(TOK_COLON); // skip the colon

        // read and store the preference list
        PreferenceList& pref_list = A[a]->get_preference_list();
        while (curtok_ != TOK_SEMICOLON) {
            std::string b = lexer_->get_lexeme();
            match(TOK_STRING);

            // if there is no vertex with that id 
            if (B.find(b) == B.end()) {
                error_occurred = true;
                std::cout << "Line " << lexer_->line_number() << ": Vertex not found: " << b << "\n";

                if (curtok_ != TOK_SEMICOLON) {
                    match(TOK_COMMA);
                }
                continue;
            }

            // if b is already present in a's pref list
            if (pref_list.find(B[b]) != pref_list.cend()) {
                error_occurred = true;
                std::cout << "Line " << lexer_->line_number() << ": Vertex " << b << " is inserted multiple times in " << a << "'s preference list\n";
            }

            pref_list.emplace_back(B[b]);

            // if there are more vertices, they must
            // be delimited using commas
            if(curtok_ == TOK_STRING) {
                std::string v = lexer_->get_lexeme();
                if(A.find(v) != A.end()) {
                    match(TOK_SEMICOLON);
                }
                else {
                    match(TOK_COMMA);
                }
            }
            else if(curtok_ == TOK_AT) {
                match(TOK_SEMICOLON);
            }
            else if (curtok_ != TOK_SEMICOLON) {
                match(TOK_COMMA);
            }
        }

        // preference list should be delimited by a semicolon
        match(TOK_SEMICOLON);
    }

    // directive should be properly terminated
    match(TOK_AT);
    match(TOK_END);
}

void GraphReader::handle_partition(BipartiteGraph::ContainerType& A, BipartiteGraph::ContainerType& B) {
    if (curtok_ == TOK_PARTITION_A) {
        match(TOK_PARTITION_A);
        read_partition(A);
    } else if(curtok_ == TOK_PARTITION_B) {
        match(TOK_PARTITION_B);
        read_partition(B);
    } else {
        std::cout << "Line " << error_message("Wrong syntax for partition", curtok_, {TOK_PARTITION_A, TOK_PARTITION_B}) << "\n";
        exit(0);
        // throw ReaderException(error_message("Wrong syntax for Partition", curtok_, {TOK_PARTITION_A, TOK_PARTITION_B}));
    }
}

void GraphReader::handle_preference_lists(BipartiteGraph::ContainerType& A, BipartiteGraph::ContainerType& B) {
    if (curtok_ == TOK_PREF_LISTS_A) {
        match(TOK_PREF_LISTS_A);
        read_preference_lists(A, B);
    } else if(curtok_ == TOK_PREF_LISTS_B) {
        match(TOK_PREF_LISTS_B);
        read_preference_lists(B, A);
    } else {
        std::cout << "Line " << error_message("Wrong syntax for preference list", curtok_, {TOK_PREF_LISTS_A, TOK_PREF_LISTS_B}) << "\n";
        exit(0);
        // throw ReaderException(error_message("Wrong syntax for Preference List", curtok_, {TOK_PREF_LISTS_A, TOK_PREF_LISTS_B}));
    }
}

void GraphReader::compatibility_checker(BipartiteGraph::ContainerType& A, BipartiteGraph::ContainerType& B)
{
    bool flag = true;
    // Checks if every 'b' in a's preference list also contains 'a' in its preference list
    for(auto it : A)
    {
        auto v = it.second;
        PreferenceList& pref_list = v->get_preference_list();
        for(auto u : pref_list)
        {
            PreferenceList& pref_list_b = (u.vertex)->get_preference_list();
            // If v is not present in pref_list_b, the input is not compatible
            if (pref_list_b.find(v) == pref_list_b.cend())
            {
                if(flag)
                {
                    std::cout << "\nIncompatible preference lists \n";
                    flag = false;
                }
                std::cout << (u.vertex)->get_id() << " doesn't have " << v->get_id() << " in its preferences but vice-versa is true." << "\n";
            }
        }
    }

    if(!flag)   exit(0);
}

std::shared_ptr<BipartiteGraph> GraphReader::read_graph() {
    BipartiteGraph::ContainerType A, B;
    
    // read the partitions
    match(TOK_AT);
    Token partition = curtok_;
    handle_partition(A, B);

    match(TOK_AT);
    // shouldn't have partition listed twice
    if (curtok_ != partition) {
        handle_partition(A, B);
    } else {
        //match((partition == TOK_PARTITION_A) ? TOK_PARTITION_B : TOK_PARTITION_A);
        std::cout << "Line " << error_message("duplicate partition listing", curtok_,
                                             {(partition == TOK_PARTITION_A) ? TOK_PARTITION_B : TOK_PARTITION_A});
        exit(0);
        // throw ReaderException(error_message("duplicate partition listing", curtok_,
        //                                     {(partition == TOK_PARTITION_A) ? TOK_PARTITION_B : TOK_PARTITION_A}));
    }

    // read the preference lists
    match(TOK_AT);
    Token pref_lists = curtok_;
    handle_preference_lists(A, B);

    match(TOK_AT);
    // shouldn't have preference lists twice
    if (curtok_ != pref_lists) {
        handle_preference_lists(A, B);
    } else {
        // match((pref_lists == TOK_PREF_LISTS_A) ? TOK_PREF_LISTS_B : TOK_PREF_LISTS_A);
        std::cout << "Line " << error_message("duplicate preference listing", curtok_,
                                             {(pref_lists == TOK_PREF_LISTS_A) ? TOK_PREF_LISTS_B : TOK_PREF_LISTS_A});
        exit(0);
        // throw ReaderException(error_message("duplicate preference listing", curtok_,
        //                                     {(pref_lists == TOK_PREF_LISTS_A) ? TOK_PREF_LISTS_B : TOK_PREF_LISTS_A}));
    }

    if(!error_occurred) {
        compatibility_checker(A, B);
        compatibility_checker(B, A);
    }

    if (error_occurred) {
        return NULL;
    }

    return std::make_shared<BipartiteGraph>(A, B);
}
