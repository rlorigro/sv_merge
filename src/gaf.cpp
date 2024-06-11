#include "gaf.hpp"

#include <unordered_map>
#include <functional>
#include <iostream>

using std::unordered_map;
using std::function;
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


void GafAlignment::parse_string_as_path(const string& p, vector<pair<string,bool> >& result){
    if (p[0] != '>' and p[0] != '<'){
        throw runtime_error("ERROR: path does not start with > or <, GAF stable path format not supported");
    }

    result.clear();

    for (auto c: p){
        if (c == '>' or c == '<'){
            result.emplace_back();
            result.back().second = parse_path_reversal_token(c);
        }
        else{
            result.back().first += c;
        }
    }
}


void GafAlignment::set_path(const string& p){
    parse_string_as_path(p, path);

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


void GafAlignment::for_step_in_path(const function<void(const string& step_name, bool is_reverse)>& f) const{
    for (const auto& [name, r]: path){
        f(name, r);
    }
}


const vector<pair<string,bool> >& GafAlignment::get_path() const{
    return path;
}


const pair<string,bool>& GafAlignment::get_step_of_path(size_t index) const{
    return path.at(index);
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


bool GafAlignment::parse_path_reversal_token(char c){
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


void GafAlignment::get_qualities(vector<uint8_t>& result){
    throw runtime_error("ERROR: get_qualities not implemented for GAF alignment");
}


void GafAlignment::get_tag_as_string(const string& name, string& result, bool allow_missing) const {
    throw runtime_error("ERROR: get_tag_as_string not implemented for GAF alignment");
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
        if (c == '\r'){
            throw runtime_error("ERROR: carriage return not supported: " + gaf_path.string());
        }

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


AlignmentSummary::AlignmentSummary(int32_t start, int32_t stop):
        start(start),
        stop(stop),
        n_match(0),
        n_mismatch(0),
        n_insert(0),
        n_delete(0)
{}


AlignmentSummary::AlignmentSummary():
        start(numeric_limits<int32_t>::max()),
        stop(numeric_limits<int32_t>::min()),
        n_match(0),
        n_mismatch(0),
        n_insert(0),
        n_delete(0)
{}


float AlignmentSummary::compute_identity() const{
    return (n_match) / (n_match + n_mismatch + n_insert + n_delete);
}


void AlignmentSummary::update(const sv_merge::CigarInterval& c, bool is_ref) {
    switch (c.code){
        case 7:
            n_match += float(c.length);    // =
            break;
        case 8:
            n_mismatch += float(c.length); // X
            break;
        case 1:
            n_insert += float(abs(c.query_stop - c.query_start));   // I
            break;
        case 2:
            n_delete += float(c.length);   // D
            break;
        default:
            break;
    }

    if (is_ref){
        auto [a,b] = c.get_forward_ref_interval();
        if (a < start){
            start = a;
        }
        if (b > stop){
            stop = b;
        }
    }
    else{
        auto [a,b] = c.get_forward_query_interval();
        if (a < start){
            start = a;
        }
        if (b > stop){
            stop = b;
        }
    }
}


GafSummary::GafSummary(const TransMap& transmap):
    transmap(transmap),
    apply_flanks(false)
{}


GafSummary::GafSummary(const VariantGraph& variant_graph, const TransMap& transmap, bool apply_flanks):
        transmap(transmap),
        apply_flanks(apply_flanks)
{
    variant_graph.graph.for_each_handle([&](const handle_t& h){
        auto l = variant_graph.graph.get_length(h);
        auto name = std::to_string(variant_graph.graph.get_id(h));
        node_lengths[name] = int32_t(l);
    });
}


void GafSummary::update_node(const string& node_name, const CigarInterval& c, bool insert){
    auto& result = ref_summaries[node_name];
    if (insert){
        result.emplace_back();
    }
    result.back().update(c, true);
}


void GafSummary::update_query(const string& query_name, const CigarInterval& c, bool insert){
    auto& result = query_summaries[query_name];
    if (insert){
        result.emplace_back();
    }
    result.back().update(c, false);
}


void compute_identity_and_coverage(const vector<AlignmentSummary>& alignments, int32_t length, pair<float,float>& result){
    AlignmentSummary total;
    int32_t bases_not_covered = 0;

    for (size_t i=0; i<alignments.size(); i++){
        const auto& a = alignments[i];

        total.n_match += a.n_match;
        total.n_mismatch += a.n_mismatch;
        total.n_insert += a.n_insert;
        total.n_delete += a.n_delete;

        if (i == 0){
            // Accumulate distance from start of sequence
            bases_not_covered += a.start;
//            cerr << "distance from start of sequence: " << a.start << '\n';
        }
        else {
            // Accumulate distance from previous alignment end
            auto distance = a.start - alignments[i-1].stop;

            if (distance < 0){
                throw runtime_error("ERROR: overlaps not resolved prior to computation of alignment identity/coverage");
            }

            bases_not_covered += distance;

//            cerr << "distance from previous alignment end: " << distance << '\n';
        }

        if (i == alignments.size() - 1){
            // Accumulate distance from end of sequence
            bases_not_covered += length - a.stop;
//            cerr << "distance from end of sequence: " << length - a.stop << '\n';
        }

//        cerr << i << ' ' << a.start << ',' << a.stop << ' ' << bases_not_covered << '\n';
//        cerr << '\n';
    }

    result.first = total.compute_identity();
    result.second = float(length - bases_not_covered) / float(length);
}


void GafSummary::for_each_ref_summary(const function<void(const string& name, int32_t length, float identity, float coverage)>& f) const{
    pair<float,float> identity_and_coverage;

    for (const auto& [name, length]: node_lengths){
        auto result = ref_summaries.find(name);

        if (result == ref_summaries.end()){
            f(name, length, 0, 0);
            continue;
        }

        const auto& alignments = result->second;

        compute_identity_and_coverage(alignments, length, identity_and_coverage);

        f(name, length, identity_and_coverage.first, identity_and_coverage.second);
    }
}


void GafSummary::for_each_query_summary(const function<void(const string& name, int32_t length, float identity, float coverage)>& f) const{
    pair<float,float> identity_and_coverage;

    int32_t length;
    transmap.for_each_read([&](const string& name, int64_t id){
        if (apply_flanks){
            coord_t c = transmap.get_flank_coord(id);
            length = c.second - c.first;
        }
        else {
            length = int32_t(transmap.get_sequence(id).size());
        }

        auto result = query_summaries.find(name);

        if (result == query_summaries.end()){
            f(name, length, 0, 0);
            return;
        }

        const auto& alignments = result->second;

        compute_identity_and_coverage(alignments, length, identity_and_coverage);

        f(name, length, identity_and_coverage.first, identity_and_coverage.second);
    });
}


HalfInterval::HalfInterval(size_t id, int32_t position, bool is_start):
        id(id),
        position(position),
        is_start(is_start)
{}


/**
 * Sweep algorithm to split intervals into intersections
 * @param alignments
 * @param is_ref
 */
void GafSummary::resolve_overlaps(vector<AlignmentSummary>& alignments) {
    if (alignments.size() <= 1){
        return;
    }

    vector <HalfInterval> half_intervals;
    vector <int32_t> positions;
    vector <unordered_set<size_t> > ids;
    vector <AlignmentSummary> result;

    for (size_t i=0; i<alignments.size(); i++){
        const auto& a = alignments[i];
        half_intervals.emplace_back(i,a.start,true);
        half_intervals.emplace_back(i,a.stop,false);
    }

    // How to sort labeled intervals with structure ((a,b), label) by start (a)
    auto left_comparator = [](const HalfInterval& a, const HalfInterval& b){
        if (a.position == b.position){
            return a.is_start > b.is_start;
        }
        else {
            return a.position < b.position;
        }
    };

    sort(half_intervals.begin(), half_intervals.end(), left_comparator);

    const auto& x = half_intervals[0];
    positions.emplace_back(x.position);
    ids.push_back({x.id});

    // With at least 2 intervals, should be guaranteed 4 half intervals in the vector
    for (size_t i=1; i<half_intervals.size(); i++){
        const auto& h = half_intervals[i];

//        cerr << "-- " << h.position << ',' << (h.is_start ? '[' : ')') << ',' << h.id << '\n';

        // Every time an element is added or removed, update the vector of intersected_intervals:
        //  1. Extend the previous interval
        //  2. Remove/add an element from the set (check if the start/stop is not identical to prev)
        if (h.is_start){
            // Only terminate the previous interval if this one does not have the same position
            if (h.position > positions.back()) {
                // Add new breakpoint
                positions.emplace_back(h.position);
                ids.emplace_back(ids.back());
            }
            // ADD item to last set
            ids.back().emplace(h.id);
        }
        else{
            // Only terminate the previous interval if this one does not have the same position
            if (h.position > positions.back()) {
                // Add new breakpoint
                positions.emplace_back(h.position);
                ids.emplace_back(ids.back());
            }
            // REMOVE item from last set
            ids.back().erase(h.id);
        }
    }

    // Finally construct new alignments that are the average of the overlapping alignments
    for (size_t i=0; i<ids.size() - 1; i++){
        if (ids[i].empty()){
            continue;
        }

        result.emplace_back();
        result.back().start = positions[i];
        result.back().stop = positions[i+1];

        auto l = float(result.back().stop - result.back().start);

        for (auto id: ids[i]){
            auto l_other = float(alignments[id].stop - alignments[id].start);

            result.back().n_match += alignments[id].n_match * (l / l_other);
            result.back().n_mismatch += alignments[id].n_mismatch * (l / l_other);
            result.back().n_insert += alignments[id].n_insert * (l / l_other);
            result.back().n_delete += alignments[id].n_delete * (l / l_other);
        }

        auto n = float(ids[i].size());
        result.back().n_match /= n;
        result.back().n_mismatch /= n;
        result.back().n_insert /= n;
        result.back().n_delete /= n;
    }

    alignments = std::move(result);
}


void GafSummary::resolve_all_overlaps() {
    for (auto& [name, alignments]: ref_summaries){
        resolve_overlaps(alignments);
    }
    for (auto& [name, alignments]: query_summaries){
        resolve_overlaps(alignments);
    }
}


void normalize_gaf_cigar(const interval_t& interval, CigarInterval& result, int32_t node_length, bool node_reversal){
    // REVERSAL:
    // path  alignment  result
    // -----------------------
    // 0     0          0
    // 0     1          1
    // 1     0          1
    // 1     1          0
    result.is_reverse = (node_reversal xor result.is_reverse);

    auto [a,b] = result.get_forward_ref_interval();

    // Compute distance from edge of interval (start of node in path)
    int32_t start = a - interval.first;
    int32_t stop = b - interval.first;

    if (not result.is_reverse){
        result.ref_start = start;
        result.ref_stop = stop;
    }
    else{
        result.ref_start = node_length - start;
        result.ref_stop = node_length - stop;
    }
}


void GafSummary::compute(const path& gaf_path){
    if (apply_flanks){
        compute_with_flanks(gaf_path);
    }
    else{
        compute_without_flanks(gaf_path);
    }
}


void GafSummary::compute_without_flanks(const path& gaf_path){
    // GAF alignments (in practice) always progress in the F orientation of the query. Reverse alignments are possible,
    // but redundant because the "reference" is actually a path, which is divided into individually reversible nodes.
    // Both minigraph and GraphAligner code do not allow for R alignments, and instead they use the path orientation.
    bool unclip_coords = false;

    for_alignment_in_gaf(gaf_path, [&](GafAlignment& alignment){
        // Skip nonsensical alignments
        if (alignment.get_map_quality() == 0){
            return;
        }

        string name;
        alignment.get_query_name(name);

        // First establish the set of intervals that defines the steps in the path (node lengths)
        // And maintain a map that will point from an interval_t to a step in the path (size_t)
        vector<interval_t> ref_intervals;
        unordered_map<interval_t,size_t> interval_to_path_index;

        int32_t x = 0;
        size_t i = 0;

        alignment.for_step_in_path([&](const string& step_name, bool is_reverse){
            auto l = int32_t(node_lengths.at(step_name));

            ref_intervals.emplace_back(x, x+l);
            interval_to_path_index.emplace(ref_intervals.back(), i);

            x += l;
            i++;
        });

        query_paths[name].emplace_back(alignment.get_path());

        // The query intervals for each alignment is not used because we do not use flanks
        vector<interval_t> query_intervals = {};

        // The iterator has a cache where they can store previously computed info for the path
        // (to avoid using unordered_map.at() for every cigar)
        string prev_node_name;
        size_t prev_path_index = -1;
        int32_t node_length = -1;
        string node_name;
        bool node_reversal;

        for_cigar_interval_in_alignment(unclip_coords, alignment, ref_intervals, query_intervals,
            // REF (node) coordinate iterator:
            [&](const CigarInterval& i, const interval_t& interval){
                // Need a copy of the cigar interval that we can normalize for stupid GAF double redundant reversal flag
                auto i_norm = i;

                auto path_index = interval_to_path_index.at(interval);

                // Don't recompute path/node stats unless we have advanced to a new node in the path
                if (path_index != prev_path_index) {
                    std::tie(node_name, node_reversal) = alignment.get_step_of_path(path_index);
                    node_length = int32_t(node_lengths.at(node_name));
                }

                normalize_gaf_cigar(interval, i_norm, node_length, node_reversal);

                // If this is a new alignment for this ref/query, inform the GafSummary to initialize a new block
                bool new_ref_alignment = node_name != prev_node_name;
                bool new_query_alignment = prev_node_name.empty();

                // Update the last alignment block in the GafSummary for ref and query
                update_node(node_name, i_norm, new_ref_alignment);
                update_query(name, i_norm, new_query_alignment);

                prev_node_name = node_name;
                prev_path_index = path_index;
            },
            // QUERY coordinate iterator: DO NOTHING HERE
            {});
    });
    resolve_all_overlaps();
}


void GafSummary::compute_with_flanks(const path& gaf_path){
    // GAF alignments (in practice) always progress in the F orientation of the query. Reverse alignments are possible,
    // but redundant because the "reference" is actually a path, which is divided into individually reversible nodes.
    // Both minigraph and GraphAligner code do not allow for R alignments, and instead they use the path orientation.
    bool unclip_coords = false;

    for_alignment_in_gaf(gaf_path, [&](GafAlignment& alignment){
        // Skip nonsensical alignments
        if (alignment.get_map_quality() == 0){
            return;
        }

        string name;
        alignment.get_query_name(name);

        // First establish the set of intervals that defines the steps in the path (node lengths)
        // And maintain a map that will point from an interval_t to a step in the path (size_t)
        vector<interval_t> ref_intervals;
        unordered_map<interval_t,size_t> interval_to_path_index;

        int32_t x = 0;
        size_t i = 0;

        alignment.for_step_in_path([&](const string& step_name, bool is_reverse){
            auto l = int32_t(node_lengths.at(step_name));

            ref_intervals.emplace_back(x, x+l);
            interval_to_path_index.emplace(ref_intervals.back(), i);

            x += l;
            i++;
        });

        query_paths[name].emplace_back(alignment.get_path());

        // The query intervals for each alignment is just the inner flank bounds
        vector<interval_t> query_intervals_i = {transmap.get_flank_coord(name)};

        // The ref and query iterators both have a cache where they can store previously computed info for the path
        // (to avoid using unordered_map.at() for every cigar)
        string prev_node_name;
        size_t prev_path_index = -1;
        int32_t node_length = -1;
        string node_name;
        bool node_reversal;
        bool new_query_alignment = true;

        for_cigar_interval_in_alignment(unclip_coords, alignment, ref_intervals, query_intervals_i,
        // REF (node) coordinate iterator:
        [&](const CigarInterval& i, const interval_t& interval){
            // Need a copy of the cigar interval that we can normalize for stupid GAF double redundant reversal flag
            auto i_norm = i;

            auto path_index = interval_to_path_index.at(interval);

            // Don't recompute path/node stats unless we have advanced to a new node in the path
            if (path_index != prev_path_index) {
                std::tie(node_name, node_reversal) = alignment.get_step_of_path(path_index);
                node_length = int32_t(node_lengths.at(node_name));
            }

            normalize_gaf_cigar(interval, i_norm, node_length, node_reversal);

            // If this is a new alignment for this ref/query, inform the GafSummary to initialize a new block
            bool new_ref_alignment = node_name != prev_node_name;

            // Update the last alignment block in the GafSummary for REF ONLY
            update_node(node_name, i_norm, new_ref_alignment);

            prev_node_name = node_name;
        },
        // QUERY coordinate iterator:
        [&](const CigarInterval& i, const interval_t& interval){
            // Update the alignment block in the GafSummary for QUERY
            // We don't care about "normalizing" the interval here because we aren't operating in ref/node space.
            // Every cigar within the query interval gets added to the summary.
            update_query(name, i, new_query_alignment);
            new_query_alignment = false;
        });
    });

    resolve_all_overlaps();
}

}
