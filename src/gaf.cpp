#include "gaf.hpp"

#include <iostream>

using std::cerr;


namespace sv_merge {


void GafAlignment::for_each_cigar_tuple(const function<void(const CigarTuple& cigar)>& f){
    for (auto& c: cigar){
        f(c);
    }
}


void GafAlignment::load_cigar(const string& cigar_string){
    string length_token;

    for (auto c: cigar_string){
        if (isalpha(c) or c == '='){
            cigar.emplace_back(stoll(length_token),cigar_char_to_code[c]);
            length_token.clear();
        }
        else{
            length_token += c;
        }
    }
}


void GafAlignment::clear_cigar(){
    cigar.clear();
}


void GafAlignment::add_tag(const string& token){
    auto i = token.find_last_of(':');

    if (i == string::npos){
        throw runtime_error("ERROR: no colon found in non-mandatory GAF column (tag)");
    }

    auto tag_type = token.substr(0,i);
    auto tag_value = token.substr(i+1,token.size() - (i+1));

    cerr << tag_type << ',' << tag_value << '\n';

    if (tag_type == "tp:A"){
        if (tag_value == "P"){
            primary = true;
        }
        else if (tag_value == "S"){
            primary = false;
        }
        else if (tag_value == "I"){
            throw runtime_error("ERROR: inversion not supported in GAF reader: " + token);
        }
        else {
            throw runtime_error("ERROR: unrecognized value in alignment tag " + token);
        }
    }

    if (tag_type == "cg:Z"){
        cigar.clear();
        load_cigar(tag_value);
    }
}


void GafAlignment::set_path(const vector<pair<string,bool> >& p){
    path = p;
}


void GafAlignment::set_path(const string& p){
    if (p[0] != '>' and p[0] != '<'){
        throw runtime_error("ERROR: path does not start with > or <, GAF stable path format not supported");
    }

    path.clear();

    for (auto c: p){
        if (c == '>' or c == '<'){
            path.emplace_back();
            path.back().second = parse_path_reversal_token(c);
        }
        else{
            path.back().first += c;
        }
    }

    // Hopefully this never happen but is technically possible according to the GAF spec
    if (reversal){
        vector<pair<string,bool> > reverse_path;
        for (auto iter=path.rbegin(); iter!=path.rend(); iter++){
            reverse_path.emplace_back(*iter);
            reverse_path.back().second = !reverse_path.back().second;
        }
        path = std::move(reverse_path);
    }
}


void GafAlignment::set_query_name(const string& name){
    query_name = name;
}


void GafAlignment::set_ref_name(const string& name){
    ref_name = name;
}


void GafAlignment::set_query_length(int32_t length){
    query_length = length;
}


void GafAlignment::set_query_start(int32_t start){
    query_start = start;
}


void GafAlignment::set_query_stop(int32_t stop){
    query_stop = stop;
}


void GafAlignment::set_path_length(int32_t length){
    path_length = length;
}


void GafAlignment::set_path_start(int32_t start){
    path_start = start;
}


void GafAlignment::set_path_stop(int32_t stop){
    path_stop = stop;
}


void GafAlignment::set_reversal(bool r){
    reversal = r;
}


void GafAlignment::set_is_primary(bool p) {
    primary = p;
}


void GafAlignment::set_n_match(int32_t n){
    n_match = n;
}


void GafAlignment::set_alignment_length(int32_t length){
    alignment_length = length;
}


void GafAlignment::set_map_quality(int32_t q){
    map_quality = q;
}


string GafAlignment::get_path_string() const{
    string p;

    for (const auto& [name, r]: path){
        p += r ? '<' : '>';
        p += name;
    }

    return p;
}


string GafAlignment::get_query_name() const{
    return query_name;
}


string GafAlignment::get_ref_name() const{
    return ref_name;
}


int32_t GafAlignment::get_query_length() const{
    return query_length;
}


int32_t GafAlignment::get_query_start() const{
    return query_start;
}


int32_t GafAlignment::get_query_stop() const{
    return query_stop;
}


int32_t GafAlignment::get_path_length() const{
    return path_length;
}


int32_t GafAlignment::get_path_start() const{
    return path_start;
}


int32_t GafAlignment::get_path_stop() const{
    return path_stop;
}


const pair<string,bool>& GafAlignment::get_path_step(int32_t index) const{
    return path[index];
}


bool GafAlignment::is_reverse() const{
    return reversal;
}


bool GafAlignment::is_primary() const{
    return primary;
}


bool GafAlignment::is_supplementary() const{
    throw runtime_error("ERROR: is_supplementary not implemented for GafAlignment");
}


int32_t GafAlignment::get_n_match() const{
    return n_match;
}


int32_t GafAlignment::get_alignment_length() const{
    return alignment_length;
}


int32_t GafAlignment::get_map_quality() const{
    return map_quality;
}


bool GafAlignment::parse_path_reversal_token(char c) const{
    bool r;

    if (c == '>'){
        r = false;
    }
    else if (c == '<'){
        r = true;
    }
    else{
        throw runtime_error("ERROR: unrecognized reversal token: " + string(1,c));
    }

    return r;
}


void GafAlignment::for_each_cigar_interval(bool unclip_coords, const function<void(const CigarInterval& cigar)>& f){
    CigarInterval c;

    // Initialize the cigar interval
    c.query_start = get_query_start();
    c.ref_start = get_ref_start();
    c.is_reverse = is_reverse();

    if (c.is_reverse){
        c.query_start = get_query_length();
    }

    for_each_cigar_tuple([&](const CigarTuple& tuple){
        c.code = tuple.code;
        c.length = tuple.length;

        // Update interval bounds for this cigar interval
        if (c.is_reverse) {
            c.query_stop = c.query_start - is_query_move[c.code]*c.length;
        }
        else{
            c.query_stop = c.query_start + is_query_move[c.code]*c.length;
        }

        c.ref_stop = c.ref_start + is_ref_move[c.code]*c.length;

        // Temporarily flip the start/stop so that it is conventionally interpretable
        c.set_query_interval_forward();

        f(c);

        // Revert to backwards intervals for iteration/update
        if (c.is_reverse){
            c.set_query_interval_reverse();
        }

        c.ref_start = c.ref_stop;
        c.query_start = c.query_stop;
    });
}


void GafAlignment::get_query_sequence(string& result){
    throw runtime_error("ERROR: get_query_sequence not implemented for GAF alignment");
}


void GafAlignment::get_query_sequence(string& result, int32_t start, int32_t stop){
    throw runtime_error("ERROR: get_query_sequence not implemented for GAF alignment");
}


void GafAlignment::get_query_name(string& result) const{
    result.clear();
    result = query_name;
}


int32_t GafAlignment::get_ref_start() const{
    return path_start;
}


int32_t GafAlignment::get_ref_stop() const{
    return path_stop;
}


bool GafAlignment::is_unmapped() const{
    // It is not possible for a PAF/GAF alignment to represent an unmapped query, so is_unmapped always returns false
    return false;
}


bool parse_reversal_token(const string& token){
    bool reversal;

    if (token == "+"){
        reversal = false;
    }
    else if (token == "-"){
        reversal = true;
    }
    else{
        throw runtime_error("ERROR: unrecognized reversal token: " + token);
    }

    return reversal;
}


// Col 	Type 	Description
// 1 	string 	Query sequence name
// 2 	int 	Query sequence length
// 3 	int 	Query start (0-based; closed)
// 4 	int 	Query end (0-based; open)
// 5 	char 	Strand relative to the path: "+" or "-"
// 6 	string 	Path matching /([><][^\s><]+(:\d+-\d+)?)+|([^\s><]+)/
// 7 	int 	Path length
// 8 	int 	Start position on the path (0-based)
// 9 	int 	End position on the path (0-based)
// 10 	int 	Number of residue matches
// 11 	int 	Alignment block length
// 12 	int 	Mapping quality (0-255; 255 for missing)
//
void for_alignment_in_gaf(const path& gaf_path, const function<void(GafAlignment& alignment)>& f){
    ifstream file(gaf_path);

    if (not (file.is_open() and file.good())){
        throw runtime_error("ERROR: could not read file: " + gaf_path.string());
    }

    char c;

    int64_t n_delimiters = 0;
    char delimiter = '\t';

    GafAlignment a;
    string token;

    while (file.get(c)){
        if (c == delimiter){
            switch (n_delimiters){
                case 0:
                    a.set_query_name(token);
                    break;
                case 1:
                    a.set_query_length(stoll(token));
                    break;
                case 2:
                    a.set_query_start(stoll(token));
                    break;
                case 3:
                    a.set_query_stop(stoll(token));
                    break;
                case 4:
                    a.set_reversal(parse_reversal_token(token));
                    break;
                case 5:
                    a.set_path(token);
                    break;
                case 6:
                    a.set_path_length(stoll(token));
                    break;
                case 7:
                    a.set_path_start(stoll(token));
                    break;
                case 8:
                    a.set_path_stop(stoll(token));
                    break;
                case 9:
                    a.set_n_match(stoll(token));
                    break;
                case 10:
                    a.set_alignment_length(stoll(token));
                    break;
                case 11:
                    a.set_map_quality(stoll(token));
                    break;
                default:
                    a.add_tag(token);
                    break;
            }

            n_delimiters++;
            token.clear();
            continue;
        }

        else if (c == '\n'){
            if (n_delimiters == 11){
                // Handle case where there is no delimiter after the last mandatory column
                a.set_map_quality(stoll(token));
            }
            if (n_delimiters > 11){
                // Handle case where there is no delimiter after the last tag
                a.add_tag(token);
            }

            f(a);
            a.clear_cigar();
            n_delimiters = 0;
            token.clear();
            continue;
        }
        else{
            if (n_delimiters == 5){
                if (c == ':') {
                    throw runtime_error("ERROR: found ':' in path, GAF stable path coordinates not supported");
                }
            }
            token += c;
        }
    }
}


void for_alignment_in_gaf(const path& gaf_path, const function<void(Alignment& alignment)>& f){
    for_alignment_in_gaf(gaf_path, [&](GafAlignment& a){
        f(a);
    });
}


}
