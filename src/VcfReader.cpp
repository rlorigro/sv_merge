#include "VcfReader.hpp"
#include "misc.hpp"

using sv_merge::equal_ignore_case;
using sv_merge::lowercase_string;

using std::min;
using std::numeric_limits;
using std::string;
using std::streamsize;

namespace sv_merge {

const char VcfReader::VCF_COMMENT = '#';
const string VcfReader::VCF_HEADER_PREFIX = "#CHROM";
const uint8_t VcfReader::VCF_HEADER_PREFIX_LENGTH = VCF_HEADER_PREFIX.length();
const uint8_t VcfReader::N_NONSAMPLE_FIELDS_VCF = 9;
const char VcfReader::LINE_END = '\n';
const char VcfReader::VCF_SEPARATOR = '\t';
const char VcfReader::VCF_MISSING_CHAR = '.';
const char VcfReader::SYMBOLIC_CHAR_OPEN = '<';
const char VcfReader::SYMBOLIC_CHAR_CLOSE = '>';
const char VcfReader::BREAKEND_CHAR_OPEN = '[';
const char VcfReader::BREAKEND_CHAR_CLOSE = ']';
const string VcfReader::PASS_STR = "PASS";
const char VcfReader::ALT_SEPARATOR = ',';
const char VcfReader::GT_SEPARATOR = ':';
const char VcfReader::BREAKEND_SEPARATOR = ':';
const char VcfReader::UNPHASED_CHAR = '/';
const char VcfReader::PHASED_CHAR = '|';
const string VcfReader::CHR_STR_LOWER = "chr";
const string VcfReader::CHR_STR_UPPER = "CHR";
const char VcfReader::UNKNOWN_BASE = 'N';

const char VcfReader::INFO_ASSIGNMENT = '=';
const char VcfReader::INFO_SEPARATOR = ';';
const string VcfReader::SVTYPE_STR = "SVTYPE";
const string VcfReader::SVLEN_STR = "SVLEN";
const string VcfReader::END_STR = "END";
const string VcfReader::CIPOS_STR = "CIPOS";
const uint16_t VcfReader::CIPOS_STR_LENGTH = CIPOS_STR.length();
const string VcfReader::CIEND_STR = "CIEND";
const uint16_t VcfReader::CIEND_STR_LENGTH = CIEND_STR.length();
const string VcfReader::CILEN_STR = "CILEN";
const uint16_t VcfReader::CILEN_STR_LENGTH = CILEN_STR.length();
const char VcfReader::CI_SEPARATOR = ',';
const string VcfReader::STDEV_POS_STR = "STDEV_POS";
const uint16_t VcfReader::STDEV_POS_STR_LENGTH = STDEV_POS_STR.length();
const string VcfReader::STDEV_LEN_STR = "STDEV_LEN";
const uint16_t VcfReader::STDEV_LEN_STR_LENGTH = STDEV_LEN_STR.length();
const string VcfReader::STDEV_END_STR = "STDEV_END";
const uint16_t VcfReader::STDEV_END_STR_LENGTH = STDEV_END_STR.length();
const string VcfReader::PRECISE_STR = "PRECISE";
const string VcfReader::IMPRECISE_STR = "IMPRECISE";
const string VcfReader::MATEID_STR = "MATEID";

const uint8_t VcfReader::TYPE_INSERTION = 1;
const uint8_t VcfReader::TYPE_DELETION = 2;
const uint8_t VcfReader::TYPE_INVERSION = 3;
const uint8_t VcfReader::TYPE_DUPLICATION = 4;
const uint8_t VcfReader::TYPE_BREAKEND = 5;
const uint8_t VcfReader::TYPE_REPLACEMENT = 6;
const uint8_t VcfReader::TYPE_CNV = 7;

const string VcfReader::DEL_STR = "DEL";
const string VcfReader::DEL_ME_STR = "DEL:ME";
const string VcfReader::INS_STR = "INS";
const string VcfReader::INS_ME_STR = "INS:ME";
const string VcfReader::INS_NOVEL_STR = "INS:NOVEL";
const string VcfReader::DUP_STR = "DUP";
const string VcfReader::DUP_TANDEM_STR = "DUP:TANDEM";
const string VcfReader::DUP_INT_STR = "DUP:INT";
const string VcfReader::INV_STR = "INV";
const string VcfReader::BND_STR = "BND";
const string VcfReader::CNV_STR = "CNV";

const uint64_t VcfRecord::STREAMSIZE_MAX = numeric_limits<streamsize>::max();


string get_vcf_name_prefix(const path& vcf){
    string name_prefix = vcf.filename().string();

    if (name_prefix.ends_with(".gz")){
        name_prefix = name_prefix.substr(0,name_prefix.size() - 3);
    }
    if (name_prefix.ends_with(".vcf")){
        name_prefix = name_prefix.substr(0,name_prefix.size() - 4);
    }

    std::replace(name_prefix.begin(), name_prefix.end(), '.', '_');

    return name_prefix;
}


VcfReader::VcfReader(const path& vcf_path, int32_t progress_n_lines, bool high_qual_only, double min_qual, bool pass_only, int32_t min_sv_length, int32_t n_samples_to_load, double min_allele_frequency, double min_nonmissing_frequency) {
    this->vcf_path=vcf_path;
    this->progress_n_lines=progress_n_lines;
    this->high_qual_only=high_qual_only;
    this->min_qual=min_qual;
    this->pass_only=pass_only;
    this->min_sv_length=min_sv_length;
    this->n_samples_to_load=n_samples_to_load;
    sample_ids=vector<string>(0);
    n_samples_in_vcf=0;
    this->min_allele_frequency=min_allele_frequency;
    this->min_nonmissing_frequency=min_nonmissing_frequency;
}


VcfReader::VcfReader(const path& vcf_path) {
    this->vcf_path=vcf_path;
    progress_n_lines=0;
    high_qual_only=false;
    min_qual=0.0;
    pass_only=false;
    min_sv_length=0;
    n_samples_to_load=INT32_MAX;
    sample_ids=vector<string>(0);
    n_samples_in_vcf=0;
    min_allele_frequency=0.0;
    min_nonmissing_frequency=0.0;
}


/**
 * @param type the value of the SVTYPE key in the INFO field;
 * @return -1 iff `type` does not represent a supported SV type.
 */
int8_t string_to_svtype_info(const string& type) {
    if (type.empty()) return -1;
    if (equal_ignore_case(type,VcfReader::DEL_STR) || equal_ignore_case(type,VcfReader::DEL_ME_STR)) return VcfReader::TYPE_DELETION;
    else if (equal_ignore_case(type,VcfReader::INS_STR) || equal_ignore_case(type,VcfReader::INS_ME_STR) || equal_ignore_case(type,VcfReader::INS_NOVEL_STR)) return VcfReader::TYPE_INSERTION;
    else if (equal_ignore_case(type,VcfReader::DUP_STR) || equal_ignore_case(type,VcfReader::DUP_TANDEM_STR) || equal_ignore_case(type,VcfReader::DUP_INT_STR)) return VcfReader::TYPE_DUPLICATION;
    else if (equal_ignore_case(type,VcfReader::INV_STR)) return VcfReader::TYPE_INVERSION;
    else if (equal_ignore_case(type,VcfReader::BND_STR)) return VcfReader::TYPE_BREAKEND;
    else if (equal_ignore_case(type,VcfReader::CNV_STR)) return VcfReader::TYPE_CNV;
    else return -1;
}


/**
 * @param type a symbolic ALT;
 * @return -1 iff `type` does not represent a supported SV type.
 */
int8_t string_to_svtype_alt(const string& type) {
    if (type.empty()) return -1;
    if (type.starts_with(VcfReader::BREAKEND_CHAR_OPEN) || type.starts_with(VcfReader::BREAKEND_CHAR_CLOSE) || type.ends_with(VcfReader::BREAKEND_CHAR_OPEN) || type.ends_with(VcfReader::BREAKEND_CHAR_CLOSE) || type==VcfReader::SYMBOLIC_CHAR_OPEN+VcfReader::BND_STR+VcfReader::SYMBOLIC_CHAR_CLOSE) return VcfReader::TYPE_BREAKEND;
    else if (equal_ignore_case(type,VcfReader::SYMBOLIC_CHAR_OPEN+VcfReader::DEL_STR+VcfReader::SYMBOLIC_CHAR_CLOSE) || equal_ignore_case(type,VcfReader::SYMBOLIC_CHAR_OPEN+VcfReader::DEL_ME_STR+VcfReader::SYMBOLIC_CHAR_CLOSE)) return VcfReader::TYPE_DELETION;
    else if (equal_ignore_case(type,VcfReader::SYMBOLIC_CHAR_OPEN+VcfReader::INS_STR+VcfReader::SYMBOLIC_CHAR_CLOSE) || equal_ignore_case(type,VcfReader::SYMBOLIC_CHAR_OPEN+VcfReader::INS_ME_STR+VcfReader::SYMBOLIC_CHAR_CLOSE) || equal_ignore_case(type,VcfReader::SYMBOLIC_CHAR_OPEN+VcfReader::INS_NOVEL_STR+VcfReader::SYMBOLIC_CHAR_CLOSE)) return VcfReader::TYPE_INSERTION;
    else if (equal_ignore_case(type,VcfReader::SYMBOLIC_CHAR_OPEN+VcfReader::DUP_STR+VcfReader::SYMBOLIC_CHAR_CLOSE) || equal_ignore_case(type,VcfReader::SYMBOLIC_CHAR_OPEN+VcfReader::DUP_TANDEM_STR+VcfReader::SYMBOLIC_CHAR_CLOSE) || equal_ignore_case(type,VcfReader::SYMBOLIC_CHAR_OPEN+VcfReader::DUP_INT_STR+VcfReader::SYMBOLIC_CHAR_CLOSE)) return VcfReader::TYPE_DUPLICATION;
    else if (equal_ignore_case(type,VcfReader::SYMBOLIC_CHAR_OPEN+VcfReader::INV_STR+VcfReader::SYMBOLIC_CHAR_CLOSE)) return VcfReader::TYPE_INVERSION;
    else if (equal_ignore_case(type,VcfReader::SYMBOLIC_CHAR_OPEN+VcfReader::CNV_STR+VcfReader::SYMBOLIC_CHAR_CLOSE)) return VcfReader::TYPE_CNV;
    else return -1;
}


VcfRecord::VcfRecord(bool high_qual_only, double min_qual, bool pass_only, int32_t min_sv_length, int32_t n_samples_to_load, double min_af, double min_nmf) {
    this->pass_only=pass_only;
    this->high_qual_only=high_qual_only;
    this->min_qual=min_qual;
    this->min_sv_length=min_sv_length;
    this->n_samples_to_load=n_samples_to_load;
    this->min_af=min_af;
    this->min_nmf=min_nmf;
    if (n_samples_to_load==INT32_MAX) {
        genotypes=vector<string>(0);
        min_n_haplotypes_alt_baseline=0.0;
        min_n_haplotypes_nonmissing_baseline=0.0;
    }
    else {
        genotypes=vector<string>(n_samples_to_load);
        min_n_haplotypes_alt_baseline=n_samples_to_load*min_af;
        min_n_haplotypes_nonmissing_baseline=n_samples_to_load*min_nmf;
    }
    tmp_buffer_1=""; tmp_buffer_2="";
    tmp_pair.first=0; tmp_pair.second=0;
}


VcfRecord VcfRecord::clone() const {
    VcfRecord out = VcfRecord(false,0.0,false,0,0,0.0,0.0);
    out.is_autosomal=this->is_autosomal;
    out.qual=this->qual;
    out.pos=this->pos;
    out.chrom=this->chrom;
    out.id=this->id;
    out.ref=this->ref;
    out.alt=this->alt;
    out.filter=this->filter;
    out.info=this->info;
    out.format=this->format;
    out.genotypes=this->genotypes;
    out.is_high_quality=this->is_high_quality;
    out.is_pass=this->is_pass;
    out.is_symbolic=this->is_symbolic;
    out.sv_type=this->sv_type;
    out.sv_length=this->sv_length;
    out.n_alts=this->n_alts;
    out.n_samples=this->n_samples;
    out.n_haplotypes_ref=this->n_haplotypes_ref;
    out.n_haplotypes_alt=this->n_haplotypes_alt;
    return out;
}


void VcfRecord::set_from_stream(ifstream& stream) {
    bool fields_skipped;
    char c;
    int32_t current_field, upper_bound;

    chrom.clear(); id.clear(); ref.clear(); alt.clear(); filter.clear(); info.clear(); format.clear();
    current_field=0; n_alts=1; n_samples=0; n_haplotypes_ref=0; n_haplotypes_alt=0; fields_skipped=false; is_autosomal=false;
    tmp_buffer_1.clear();
    while (stream.get(c)) {
        if (c==VcfReader::VCF_SEPARATOR || c==VcfReader::LINE_END) {
            fields_skipped=set_field(tmp_buffer_1,current_field,high_qual_only,min_qual,pass_only,min_sv_length,min_af,min_nmf,stream,tmp_buffer_2,tmp_pair);
            if (fields_skipped) break;
            else {
                tmp_buffer_1.clear(); current_field++;
                if (c==VcfReader::LINE_END) break;
                upper_bound=n_samples_to_load-n_samples;
                if (is_autosomal || chrom=="chrX" || chrom=="X") upper_bound<<=1;
                if ( n_haplotypes_alt+upper_bound<min_n_haplotypes_alt ||
                     (n_haplotypes_ref+n_haplotypes_alt)+upper_bound<min_n_haplotypes_nonmissing
                     ) {
                    stream.ignore(STREAMSIZE_MAX,VcfReader::LINE_END);
                    fields_skipped=true;
                    break;
                }
            }
        }
        else {
            tmp_buffer_1+=c;
            if (current_field==4 && c==VcfReader::ALT_SEPARATOR) n_alts++;
        }
    }
    if (!fields_skipped && !tmp_buffer_1.empty()) set_field(tmp_buffer_1,current_field,high_qual_only,min_qual,pass_only,min_sv_length,min_af,min_nmf,stream,tmp_buffer_2,tmp_pair);
}


bool VcfRecord::set_field(const string& field, int32_t field_id, bool high_qual_only, double min_qual, bool pass_only, int32_t min_sv_length, double min_af, double min_nmf, ifstream& stream, string& tmp_buffer, pair<uint8_t, uint8_t>& tmp_pair) {
    if (field_id==0) {
        chrom=field;
        is_autosomal = chrom!="chrX" && chrom!="X" &&  chrom!="chrY" && chrom!="Y" && chrom!="chrM" && chrom!="M";
        min_n_haplotypes_alt=ceil((is_autosomal?2.0:1.0)*min_n_haplotypes_alt_baseline);
        min_n_haplotypes_nonmissing=ceil((is_autosomal?2.0:1.0)*min_n_haplotypes_nonmissing_baseline);
    }
    else if (field_id==1) pos=stoi(field);
    else if (field_id==2) id.append(field);
    else if (field_id==3) ref.append(field);
    else if (field_id==4) {
        alt.append(field);
        is_symbolic=is_alt_symbolic();
        if (n_alts>1) {
            stream.ignore(STREAMSIZE_MAX,VcfReader::LINE_END);
            return true;
        }
        if ( equal_ignore_case(field,VcfReader::SYMBOLIC_CHAR_OPEN+VcfReader::INS_STR+VcfReader::SYMBOLIC_CHAR_CLOSE) ||
             equal_ignore_case(field,VcfReader::SYMBOLIC_CHAR_OPEN+VcfReader::INS_ME_STR+VcfReader::SYMBOLIC_CHAR_CLOSE) ||
             equal_ignore_case(field,VcfReader::SYMBOLIC_CHAR_OPEN+VcfReader::INS_NOVEL_STR+VcfReader::SYMBOLIC_CHAR_CLOSE)
             ) {
            sv_type=VcfReader::TYPE_INSERTION;
            stream.ignore(STREAMSIZE_MAX,VcfReader::LINE_END);
            return true;
        }
    }
    else if (field_id==5) {
        if (field.starts_with(VcfReader::VCF_MISSING_CHAR)) { qual=-1; is_high_quality=false; }
        else {
            qual=stof(field);
            is_high_quality=qual>=min_qual;
            if (high_qual_only && !is_high_quality) {
                stream.ignore(STREAMSIZE_MAX,VcfReader::LINE_END);
                return true;
            }
        }
    }
    else if (field_id==6) {
        filter.append(field);
        if (field.starts_with(VcfReader::VCF_MISSING_CHAR)) is_pass=false;
        else {
            is_pass=filter==VcfReader::PASS_STR;
            if (pass_only && !is_pass) {
                stream.ignore(STREAMSIZE_MAX,VcfReader::LINE_END);
                return true;
            }
        }
    }
    else if (field_id==7) {
        info.append(field);
        set_sv_type(tmp_buffer);
        set_sv_length(tmp_buffer);
        const auto alt_length = (int32_t)alt.length();
        if (sv_type==VcfReader::TYPE_INSERTION && !is_symbolic && sv_length!=alt_length-1) sv_length=alt_length-1;
        if (sv_type==-1 || sv_length<min_sv_length || sv_length==0) {
            stream.ignore(STREAMSIZE_MAX,VcfReader::LINE_END);
            return true;
        }
    }
    else if (field_id==8) {
        format.append(field);
        if (n_samples_to_load==0) {
            stream.ignore(STREAMSIZE_MAX,VcfReader::LINE_END);
            return true;
        }
    }
    else {
        n_samples++;
        if (n_samples>n_samples_to_load) {
            stream.ignore(STREAMSIZE_MAX,VcfReader::LINE_END);
            return true;
        }
        else {
            genotypes[n_samples-1].clear();
            genotypes[n_samples-1].append(field);
            ncalls_in_sample(field,tmp_pair);
            n_haplotypes_ref+=tmp_pair.first; n_haplotypes_alt+=tmp_pair.second;
        }
    }
    return false;
}


void VcfRecord::set_sv_type(string& tmp_buffer) {
    tmp_buffer.clear();
    const bool found = get_info_field(VcfReader::SVTYPE_STR,0,tmp_buffer)!=string::npos;
    if (found) sv_type=string_to_svtype_info(tmp_buffer);
    else if (alt.starts_with(VcfReader::SYMBOLIC_CHAR_OPEN) && alt.ends_with(VcfReader::SYMBOLIC_CHAR_CLOSE)) sv_type=string_to_svtype_alt(alt);
    else if (alt.starts_with(VcfReader::BREAKEND_CHAR_OPEN) || alt.starts_with(VcfReader::BREAKEND_CHAR_CLOSE) || alt.ends_with(VcfReader::BREAKEND_CHAR_OPEN) || alt.ends_with(VcfReader::BREAKEND_CHAR_CLOSE)) sv_type=string_to_svtype_alt(alt);
    else if (!alt.starts_with(VcfReader::VCF_MISSING_CHAR)) {
        const size_t ref_length = ref.length();
        const size_t alt_length = alt.length();
        if (ref_length==1 && alt_length>ref_length) sv_type=VcfReader::TYPE_INSERTION;
        else if (alt_length==1 && ref_length>alt_length) sv_type=VcfReader::TYPE_DELETION;
        else if (ref_length>1 && alt_length>1) sv_type=VcfReader::TYPE_REPLACEMENT;
        else sv_type=-1;
    }
    else sv_type=-1;
}


void VcfRecord::set_sv_length(string& tmp_buffer) {
    tmp_buffer.clear();
    if (get_info_field(VcfReader::SVLEN_STR,0,tmp_buffer)!=string::npos) {
        const int32_t length = stoi(tmp_buffer);
        sv_length=length<0?-length:length;
        if (sv_length>0) return;
    }
    if (get_info_field(VcfReader::END_STR,0,tmp_buffer)!=string::npos) {
        const int32_t end = stoi(tmp_buffer);
        sv_length=end-pos;
        if (sv_length>0) return;
    }
    if (alt.starts_with(VcfReader::SYMBOLIC_CHAR_OPEN) && alt.ends_with(VcfReader::SYMBOLIC_CHAR_CLOSE)) sv_length=INT32_MAX;
    else if (alt.starts_with(VcfReader::BREAKEND_CHAR_OPEN) || alt.starts_with(VcfReader::BREAKEND_CHAR_CLOSE) || alt.ends_with(VcfReader::BREAKEND_CHAR_OPEN) || alt.ends_with(VcfReader::BREAKEND_CHAR_CLOSE)) sv_length=INT32_MAX;
    else {
        const size_t ref_length = ref.length();
        const size_t alt_length = alt.length();
        if (ref_length==1 && alt_length>ref_length) sv_length=(int32_t)(alt_length-1);  // INS
        else if (alt_length==1 && ref_length>alt_length) sv_length=(int32_t)(ref_length-1);  // DEL
        else sv_length=(int32_t)(ref_length-1);  // Replacement
    }
}


void VcfRecord::print(ostream& stream) const {
    stream << chrom << VcfReader::VCF_SEPARATOR << pos << VcfReader::VCF_SEPARATOR << id << VcfReader::VCF_SEPARATOR << ref << VcfReader::VCF_SEPARATOR << alt << VcfReader::VCF_SEPARATOR;
    if (qual==-1) stream << '.';
    else stream << qual;
    stream << VcfReader::VCF_SEPARATOR << filter << VcfReader::VCF_SEPARATOR << info << VcfReader::VCF_SEPARATOR << format;

    for (const auto& item: genotypes) {
        stream << VcfReader::VCF_SEPARATOR << item;
    }
}


uint8_t VcfRecord::get_gt(int32_t sample, pair<int8_t,int8_t>& out) {
    if (sample>=n_samples) return 0;
    char c;
    const size_t LENGTH = genotypes[sample].length();
    size_t i;

    out.first=-1; out.second=-1;
    tmp_buffer_1.clear();
    for (i=0; i<LENGTH; i++) {
        c=genotypes[sample].at(i);
        if (c==VcfReader::PHASED_CHAR || c==VcfReader::UNPHASED_CHAR) {
            if (!tmp_buffer_1.starts_with(VcfReader::VCF_MISSING_CHAR)) out.first=(int8_t)stoi(tmp_buffer_1);
            tmp_buffer_1.clear();
        }
        else if (c==VcfReader::GT_SEPARATOR) {
            if (!tmp_buffer_1.starts_with(VcfReader::VCF_MISSING_CHAR)) out.second=(int8_t)stoi(tmp_buffer_1);
            return 2;
        }
        else tmp_buffer_1+=c;
    }
    if (!tmp_buffer_1.starts_with(VcfReader::VCF_MISSING_CHAR)) out.second=(int8_t)stoi(tmp_buffer_1);
    if (out.first==-1 && out.second!=-1) { out.first=out.second; out.second=-1; return 1; }
    else return 2;
}


bool VcfRecord::passes_constraints() const {
    return n_alts==1 &&
           (is_high_quality || !high_qual_only) &&
           (!pass_only || is_pass || filter.starts_with(VcfReader::VCF_MISSING_CHAR)) &&
           sv_type!=-1 && sv_length>=min_sv_length &&
           (n_samples_to_load==0 || (n_samples==n_samples_to_load && n_haplotypes_alt>=min_n_haplotypes_alt && n_haplotypes_ref+n_haplotypes_alt>=min_n_haplotypes_nonmissing));
}


size_t VcfRecord::get_info_field(const string &key, const size_t from, string &out) const {
    const uint8_t KEY_LENGTH = key.length()+1;
    size_t i, p, q;

    p=info.find(key,from);
    if (p==string::npos) {
        string lower = key;
        lowercase_string(lower);
        p=info.find(lower,from);
    }
    if (p==string::npos) return string::npos;
    if (key==VcfReader::END_STR && p>=2 && toupper(info.at(p-2))=='C' && toupper(info.at(p-1))=='I') {
        p=info.find(key,p+KEY_LENGTH);
        if (p==string::npos) return string::npos;
    }
    q=info.find(VcfReader::INFO_SEPARATOR,p+KEY_LENGTH);
    if (q==string::npos) q=info.length();
    out.clear();
    for (i=p+KEY_LENGTH; i<q; i++) out+=info.at(i);
    return p;
}


void VcfRecord::ncalls_in_sample(const string& buffer, pair<uint8_t, uint8_t>& out) {
    char c;
    const size_t LENGTH = buffer.length();
    size_t i, j;

    out.first=0; out.second=0;
    j=0;
    for (i=0; i<LENGTH; i++) {
        c=buffer.at(i);
        if (c==VcfReader::PHASED_CHAR || c==VcfReader::UNPHASED_CHAR) {
            if (buffer.starts_with('0')) out.first=1;
            else if (!buffer.starts_with(VcfReader::VCF_MISSING_CHAR)) out.second=1;
            j=i+1;
        }
        else if (c==VcfReader::GT_SEPARATOR) {
            c=buffer.at(j);
            if (c=='0') out.first++;
            else if (c!=VcfReader::VCF_MISSING_CHAR) out.second++;
            return;
        }
    }
    c=buffer.at(j);
    if (c=='0') out.first++;
    else if (c!=VcfReader::VCF_MISSING_CHAR) out.second++;
}


bool VcfRecord::is_phased(const string& buffer) {
    const size_t LENGTH = buffer.length();

    for (size_t i=0; i<LENGTH; i++) {
        switch (buffer.at(i)) {
            case VcfReader::PHASED_CHAR: return true;
            case VcfReader::UNPHASED_CHAR: return false;
            case VcfReader::GT_SEPARATOR: break;
        }
    }
    return true;
}


void VcfRecord::info2map(unordered_map<string,string>& map) {
    bool found;
    char c;
    uint16_t i;
    const uint16_t LENGTH = info.length();

    map.clear(); found=false;
    tmp_buffer_1.clear(); tmp_buffer_2.clear();
    for (i=0; i<LENGTH; i++) {
        c=info.at(i);
        if (c==VcfReader::INFO_ASSIGNMENT) found=true;
        else if (c==VcfReader::INFO_SEPARATOR) {
            map.emplace(tmp_buffer_1,tmp_buffer_2);
            tmp_buffer_1.clear(); tmp_buffer_2.clear();
            found=false;
        }
        else if (!found) tmp_buffer_1+=c;
        else tmp_buffer_2+=c;
    }
    if (found) map.emplace(tmp_buffer_1,tmp_buffer_2);
}


void VcfRecord::get_samples_with_alt(unordered_set<int32_t>& out) {
    int32_t i;
    const uint32_t SIZE = genotypes.size();

    out.clear();
    for (i=0; i<SIZE; i++) {
        tmp_buffer_1.clear(); tmp_buffer_1.append(genotypes.at(i));
        ncalls_in_sample(tmp_buffer_1,tmp_pair);
        if (tmp_pair.second!=0) out.emplace(i);
    }
}


void VcfRecord::get_samples_with_alt(set<int32_t>& out) {
    int32_t i;
    const uint32_t SIZE = genotypes.size();

    out.clear();
    for (i=0; i<SIZE; i++) {
        tmp_buffer_1.clear(); tmp_buffer_1.append(genotypes.at(i));
        ncalls_in_sample(tmp_buffer_1,tmp_pair);
        if (tmp_pair.second!=0) out.emplace(i);
    }
}


bool VcfRecord::is_alt_symbolic() const { return alt.at(0)==VcfReader::SYMBOLIC_CHAR_OPEN && alt.at(alt.length()-1)==VcfReader::SYMBOLIC_CHAR_CLOSE; }


uint8_t VcfRecord::is_breakend_single() const {
    const char first = alt.at(0);
    const char last = alt.at(alt.length()-1);
    if ((first==VcfReader::VCF_MISSING_CHAR&&last!=VcfReader::VCF_MISSING_CHAR) || (first!=VcfReader::VCF_MISSING_CHAR&&last==VcfReader::VCF_MISSING_CHAR)) {
        const uint32_t length = alt.length();
        if (length==2) return 0;
        else if (length>2) return 1;
    }
    return 2;
}


bool VcfRecord::is_breakend_virtual(const unordered_map<string,string>& chromosomes) {
    if (is_symbolic) return false;
    if (pos==0 || pos==chromosomes.at(chrom).length()+1) return true;
    const int32_t p = get_breakend_pos();
    if (p==INT32_MAX) return false;
    string chr_trans;
    get_breakend_chromosome(chr_trans);
    return (p==0 || p==chromosomes.at(chr_trans).length()+1);
}


void VcfRecord::get_breakend_chromosome(string& out) const {
    if (is_symbolic) { out.clear(); return; }
    char c;
    uint16_t i, p;
    const uint16_t LENGTH = alt.length();

    out.clear(); p=UINT16_MAX;
    for (i=0; i<LENGTH; i++) {
        c=alt.at(i);
        if (p==UINT16_MAX && (c==VcfReader::BREAKEND_CHAR_OPEN || c==VcfReader::BREAKEND_CHAR_CLOSE)) p=i;
        else if (c==VcfReader::BREAKEND_SEPARATOR) break;
        else if (p!=UINT16_MAX) out.push_back(c);
    }
    if (out.starts_with(VcfReader::CHR_STR_UPPER)) out.replace(0,VcfReader::CHR_STR_LOWER.length(),VcfReader::CHR_STR_LOWER);
}


int32_t VcfRecord::get_breakend_pos() {
    if (is_symbolic) return INT32_MAX;
    char c;
    uint16_t i, p;
    const uint16_t LENGTH = alt.length();

    tmp_buffer_1.clear(); p=UINT16_MAX;
    for (i=0; i<LENGTH; i++) {
        c=alt.at(i);
        if (p==UINT16_MAX) {
            if (c==VcfReader::BREAKEND_SEPARATOR) p=i;
        }
        else {
            if (c==VcfReader::BREAKEND_CHAR_OPEN || c==VcfReader::BREAKEND_CHAR_CLOSE) break;
            tmp_buffer_1.push_back(c);
        }
    }
    return tmp_buffer_1.empty()?INT32_MAX:stoi(tmp_buffer_1);
}


uint8_t VcfRecord::get_breakend_orientation_cis() const {
    if (is_symbolic) return 0;
    const char c = alt.at(0);
    return (c!=VcfReader::BREAKEND_CHAR_OPEN && c!=VcfReader::BREAKEND_CHAR_CLOSE && c!=VcfReader::VCF_MISSING_CHAR)?1:2;
}


uint8_t VcfRecord::get_breakend_orientation_trans() const {
    if (is_symbolic || is_breakend_single()<2) return 0;
    char c;

    c=alt.at(0);
    if (c==VcfReader::BREAKEND_CHAR_OPEN) return 2;
    else if (c==VcfReader::BREAKEND_CHAR_CLOSE) return 1;
    c=alt.at(alt.length()-1);
    if (c==VcfReader::BREAKEND_CHAR_OPEN) return 2;
    else if (c==VcfReader::BREAKEND_CHAR_CLOSE) return 1;
    return 0;
}


void VcfRecord::get_breakend_inserted_sequence(string& out) const {
    if (is_symbolic) { out.clear(); return; }
    const uint16_t LENGTH = alt.length();
    char c;
    uint16_t p;

    out.clear(); p=UINT16_MAX;
    for (uint16_t i=0; i<LENGTH; i++) {
        c=alt.at(i);
        if (c!=VcfReader::BREAKEND_CHAR_OPEN && c!=VcfReader::BREAKEND_CHAR_CLOSE) continue;
        if (p==UINT16_MAX) {
            if (i>1) { out.append(alt.substr(1,i-1)); return; }
            p=i;
        }
        else if (i<LENGTH-2) { out.append(alt.substr(i+1,LENGTH-i-2)); return; }
    }
    if (LENGTH>2 && (alt.at(0)==VcfReader::VCF_MISSING_CHAR || alt.at(LENGTH-1)==VcfReader::VCF_MISSING_CHAR)) { out.append(alt.substr(1,LENGTH-2)); return; }
}


void VcfRecord::get_alt_of_breakend_mate(const unordered_map<string,string>& chromosomes, string& out) {
    const int32_t pos2 = get_breakend_pos();
    get_breakend_chromosome(tmp_buffer_1);
    const char char2 = chromosomes.at(tmp_buffer_1).at(pos2-1);

    out.clear();
    const char c = alt.at(0);
    if (c==VcfReader::BREAKEND_CHAR_OPEN) {
        out.push_back(VcfReader::BREAKEND_CHAR_OPEN);
        out.append(chrom);
        out.push_back(VcfReader::BREAKEND_SEPARATOR);
        out.append(std::to_string(pos));
        out.push_back(VcfReader::BREAKEND_CHAR_OPEN);
        out.push_back(char2);
    }
    else if (c==VcfReader::BREAKEND_CHAR_CLOSE) {
        out.push_back(char2);
        out.push_back(VcfReader::BREAKEND_CHAR_OPEN);
        out.append(chrom);
        out.push_back(VcfReader::BREAKEND_SEPARATOR);
        out.append(std::to_string(pos));
        out.push_back(VcfReader::BREAKEND_CHAR_OPEN);
    }
    else if (alt.find(VcfReader::BREAKEND_CHAR_OPEN)!=string::npos) {
        out.push_back(VcfReader::BREAKEND_CHAR_CLOSE);
        out.append(chrom);
        out.push_back(VcfReader::BREAKEND_SEPARATOR);
        out.append(std::to_string(pos));
        out.push_back(VcfReader::BREAKEND_CHAR_CLOSE);
        out.push_back(char2);
    }
    else {
        out.push_back(char2);
        out.push_back(VcfReader::BREAKEND_CHAR_CLOSE);
        out.append(chrom);
        out.push_back(VcfReader::BREAKEND_SEPARATOR);
        out.append(std::to_string(pos));
        out.push_back(VcfReader::BREAKEND_CHAR_CLOSE);
    }
}


bool VcfRecord::is_precise() {
    if (get_info_field(VcfReader::IMPRECISE_STR,0,tmp_buffer_1)!=string::npos) return false;
    for (uint8_t i=0; i<=2; i++) {
        get_confidence_interval(i,tmp_pair_2);
        if (tmp_pair_2.first!=0.0 || tmp_pair_2.second!=0.0) return false;
    }
    return true;
}


void VcfRecord::get_confidence_interval(uint8_t which, pair<double,double>& out) {
    const uint8_t SIGMA_MULTIPLE = 3;  // 2 or 3 captures most of a Gaussian
    string KEY_1, KEY_2;
    char c;
    double sigma;
    uint16_t length;
    size_t i, p;

    switch (which) {
        case 0: KEY_1=VcfReader::CIPOS_STR; KEY_2=VcfReader::STDEV_POS_STR; break;
        case 1: KEY_1=VcfReader::CILEN_STR; KEY_2=VcfReader::STDEV_LEN_STR; break;
        case 2: KEY_1=VcfReader::CIEND_STR; KEY_2=VcfReader::STDEV_END_STR; break;
    }
    out.first=0.0; out.second=0.0;
    p=get_info_field(KEY_1,0,tmp_buffer_1);
    if (p!=string::npos) {
        tmp_buffer_2.clear(); length=tmp_buffer_1.length();
        for (i=0; i<length; i++) {
            c=tmp_buffer_1.at(i);
            if (c==VcfReader::CI_SEPARATOR) {
                out.first=stof(tmp_buffer_2);
                tmp_buffer_2.clear();
            }
            else tmp_buffer_2+=c;
        }
        out.second=stof(tmp_buffer_2);
        return;
    }
    p=get_info_field(KEY_2,0,tmp_buffer_1);
    if (p!=string::npos) {
        sigma=stof(tmp_buffer_1)*SIGMA_MULTIPLE;
        out.first=-sigma; out.second=sigma;
    }
}


void VcfRecord::get_confidence_interval_pos(pair<double,double>& out) { get_confidence_interval(0,out); }

void VcfRecord::get_confidence_interval_length(pair<double,double>& out) { get_confidence_interval(1,out); }

void VcfRecord::get_confidence_interval_end(pair<double,double>& out) { get_confidence_interval(2,out); }


void VcfRecord::get_reference_coordinates(bool use_confidence_intervals, coord_t& out) {
    if (sv_type==VcfReader::TYPE_INSERTION) {
        if (use_confidence_intervals) {
            get_confidence_interval_pos(tmp_pair_2);
            out.first=floor(pos+tmp_pair_2.first);
            out.second=ceil(pos+tmp_pair_2.second);
        }
        else { out.first=pos; out.second=out.first; }
    }
    else if (sv_type==VcfReader::TYPE_BREAKEND) {
        if (pos==0) {  // Virtual telomeric breakend
            out.first=INT32_MAX; out.second=INT32_MAX;
        }
        else if (use_confidence_intervals) {
            get_confidence_interval_pos(tmp_pair_2);
            out.first=floor(pos-1+tmp_pair_2.first);
            out.second=ceil(pos-1+tmp_pair_2.second);
        }
        else { out.first=pos-1; out.second=out.first; }
    }
    else if (sv_type==VcfReader::TYPE_DELETION || sv_type==VcfReader::TYPE_INVERSION || sv_type==VcfReader::TYPE_DUPLICATION || sv_type==VcfReader::TYPE_CNV) {
        if (use_confidence_intervals) {
            get_confidence_interval_pos(tmp_pair_2);
            out.first=floor(pos+tmp_pair_2.first);
            if (sv_length==INT32_MAX) out.second=INT32_MAX;
            else {
                get_confidence_interval_length(tmp_pair_2);
                out.second=ceil(pos+sv_length+tmp_pair_2.second);
            }
        }
        else {
            out.first=pos;
            out.second=sv_length==INT32_MAX?INT32_MAX:pos+sv_length;
        }
    }
    else if (sv_type==VcfReader::TYPE_REPLACEMENT) {
        out.first=pos;
        out.second=pos+sv_length;
    }
    else { out.first=INT32_MAX; out.second=INT32_MAX; }
}


bool VcfRecord::is_valid_bnd_alt() {
    if (is_symbolic) return true;
    const bool is_single = is_breakend_single();
    bool is_nonsingle;
    size_t p;

    p=alt.find(VcfReader::BREAKEND_CHAR_OPEN);
    if (p!=string::npos) is_nonsingle=alt.find(VcfReader::BREAKEND_CHAR_OPEN,p+1)!=string::npos;
    else {
        p=alt.find(VcfReader::BREAKEND_CHAR_CLOSE);
        if (p!=string::npos) is_nonsingle=alt.find(VcfReader::BREAKEND_CHAR_CLOSE,p+1)!=string::npos;
        else is_nonsingle=false;
    }
    return is_single!=is_nonsingle;
}


void VcfReader::for_record_in_vcf(const function<void(VcfRecord& record)>& callback) {
    bool is_header;
    char c;
    uint8_t i;
    int32_t n_fields;
    int64_t n_lines;
    VcfRecord record(high_qual_only,min_qual,pass_only,min_sv_length,n_samples_to_load,min_allele_frequency,min_nonmissing_frequency);
    ifstream file;
    string buffer;

    file.open(vcf_path,std::ifstream::in);
    if (!file.is_open() || !file.good()) throw runtime_error("ERROR: could not read file: " + vcf_path);
    n_lines=0; vcf_header.clear();
    while (!file.eof()) {
        c=(char)file.peek();
        if (c==VCF_COMMENT) {
            is_header=true;
            for (i=0; i<VCF_HEADER_PREFIX_LENGTH; i++) {
                file.get(c);
                if (toupper(c)!=VCF_HEADER_PREFIX.at(i)) { is_header=false; break; }
                else vcf_header.push_back(c);
            }
            if (!is_header) {
                vcf_header.clear();
                file.ignore(VcfRecord::STREAMSIZE_MAX,LINE_END);
                continue;
            }
            n_fields=0; buffer.clear();
            while (file.get(c)) {
                if (c==LINE_END || c==EOF) {
                    n_fields++;
                    if (n_fields>N_NONSAMPLE_FIELDS_VCF) sample_ids.push_back(buffer);
                    buffer.clear();
                    break;
                }
                else if (c==VCF_SEPARATOR) {
                    n_fields++;
                    if (n_fields>N_NONSAMPLE_FIELDS_VCF) sample_ids.push_back(buffer);
                    buffer.clear();
                    vcf_header.push_back(c);
                }
                else {
                    buffer+=c;
                    vcf_header.push_back(c);
                }
            }
            n_samples_in_vcf=n_fields-N_NONSAMPLE_FIELDS_VCF;
            n_samples_to_load=min(n_samples_to_load,n_samples_in_vcf);
            record.min_n_haplotypes_alt_baseline=n_samples_to_load*min_allele_frequency;
            record.min_n_haplotypes_nonmissing_baseline=n_samples_to_load*min_nonmissing_frequency;
            record.n_samples_to_load=n_samples_to_load;
            if (n_samples_to_load>record.genotypes.size()) record.genotypes.resize(n_samples_to_load);
            continue;
        }
        else if (c==EOF) break;
        record.set_from_stream(file);
        if (record.passes_constraints()) callback(record);
        n_lines++;
        if (progress_n_lines!=0 && n_lines%progress_n_lines==0) cerr << "Scanned " << n_lines << " records\n";
    }
    file.close();
    if (progress_n_lines!=0) cerr << "Scanned " << n_lines << " total records\n";
}


void VcfReader::print_minimal_header(ofstream& out) const {
    out << VCF_FORMAT_LINE << '\n';
    out << vcf_header << '\n';
}


void VcfReader::get_minimal_header(string& out) const {
    out=VCF_FORMAT_LINE+'\n'+vcf_header;
}


void VcfReader::get_file_path(path& out) const {
    out=vcf_path;
}

}