#include "VcfReader.hpp"

using std::string;
using std::numeric_limits;
using std::streamsize;

namespace sv_merge {


/**
 * @param type the value of the SVTYPE key in the INFO field;
 * @return -1 iff `type` does not represent a supported SV type.
 */
int8_t string_to_svtype_info(const string& type) {
    if (type.length()==0) return -1;
    if (type==DEL_STR || type==DEL_ME_STR) return TYPE_DELETION;
    else if (type==INS_STR || type==INS_ME_STR || type==INS_NOVEL_STR) return TYPE_INSERTION;
    else if (type==DUP_STR || type==DUP_TANDEM_STR || type==DUP_INT_STR) return TYPE_DUPLICATION;
    else if (type==INV_STR) return TYPE_INVERSION;
    else if (type==BND_STR) return TYPE_BREAKEND;
    else return -1;
}


/**
 * @param type a symbolic ALT;
 * @return -1 iff `type` does not represent a supported SV type.
 */
int8_t string_to_svtype_alt(const string& type) {
    if (type.length()==0) return -1;
    if (type.starts_with('[') || type.starts_with(']') || type.ends_with('[') || type.ends_with(']') || type=="<"+BND_STR+">") return TYPE_BREAKEND;
    else if (type=="<"+DEL_STR+">" || type=="<"+DEL_ME_STR+">") return TYPE_DELETION;
    else if (type=="<"+INS_STR+">" || type=="<"+INS_ME_STR+">" || type=="<"+INS_NOVEL_STR+">") return TYPE_INSERTION;
    else if (type=="<"+DUP_STR+">" || type=="<"+DUP_TANDEM_STR+">" || type=="<"+DUP_INT_STR+">") return TYPE_DUPLICATION;
    else if (type=="<"+INV_STR+">") return TYPE_INVERSION;
    else return -1;
}


VcfRecord::VcfRecord(const string& chrom, bool is_diploid, bool high_qual_only, float min_qual, bool pass_only, uint32_t min_sv_length, bool skip_samples, uint32_t n_samples_in_vcf, float min_allele_frequency, float min_nonmissing_frequency) {
    PASS_ONLY=pass_only;
    HIGH_QUAL_ONLY=high_qual_only;
    MIN_QUAL=min_qual;
    SKIP_SAMPLES=skip_samples;
    MIN_SV_LENGTH=min_sv_length;
    this->chrom=chrom;
    this->is_diploid=is_diploid;
    N_SAMPLES_IN_VCF=n_samples_in_vcf;
    genotypes=vector<string>(N_SAMPLES_IN_VCF);
    MIN_N_HAPLOTYPES_ALT=ceil((is_diploid?2.0:1.0)*N_SAMPLES_IN_VCF*min_allele_frequency);
    MIN_N_HAPLOTYPES_NONMISSING=ceil((is_diploid?2.0:1.0)*N_SAMPLES_IN_VCF*min_nonmissing_frequency);
}


void VcfRecord::set(ifstream& stream, string& tmp_buffer_1, string& tmp_buffer_2, pair<uint8_t, uint8_t>& tmp_pair) {
    bool fields_skipped;
    char c;
    uint32_t current_field, upper_bound;

    id.clear(); ref.clear(); alt.clear(); filter.clear(); info.clear(); format.clear(); genotypes.clear();
    stream.ignore(numeric_limits<streamsize>::max(),VCF_SEPARATOR);  // Skipping CHROM
    current_field=1; n_alts=1; n_samples=0; n_haplotypes_ref=0; n_haplotypes_alt=0; fields_skipped=false;
    tmp_buffer_1.clear();
    while (stream.get(c)) {
        if (c==VCF_SEPARATOR || c==LINE_END) {
            fields_skipped=set_field(tmp_buffer_1,current_field,HIGH_QUAL_ONLY,MIN_QUAL,PASS_ONLY,MIN_SV_LENGTH,SKIP_SAMPLES,stream,tmp_buffer_2,tmp_pair);
            if (fields_skipped) break;
            else {
                tmp_buffer_1.clear(); current_field++;
                if (c==LINE_END) break;
                upper_bound=N_SAMPLES_IN_VCF-n_samples;
                if (is_diploid) upper_bound<<=1;
                if ( n_haplotypes_alt+upper_bound<MIN_N_HAPLOTYPES_ALT ||
                     (n_haplotypes_ref+n_haplotypes_alt)+upper_bound<MIN_N_HAPLOTYPES_NONMISSING
                     ) {
                    stream.ignore(numeric_limits<streamsize>::max(),LINE_END);
                    fields_skipped=true;
                    break;
                }
            }
        }
        else {
            tmp_buffer_1+=c;
            if (current_field==4 && c==ALT_SEPARATOR) n_alts++;
        }
    }
    if (!fields_skipped && tmp_buffer_1.length()!=0) set_field(tmp_buffer_1,current_field,HIGH_QUAL_ONLY,MIN_QUAL,PASS_ONLY,MIN_SV_LENGTH,SKIP_SAMPLES,stream,tmp_buffer_2,tmp_pair);
}


bool VcfRecord::set_field(const string& field, uint32_t field_id, bool high_qual_only, float min_qual, bool pass_only, uint32_t min_sv_length, bool skip_samples, ifstream& stream, string& tmp_buffer, pair<uint8_t, uint8_t>& tmp_pair) {
    // field_id is never zero since CHROM is skipped
    if (field_id==1) pos=stoul(field);
    else if (field_id==2) id+=field;
    else if (field_id==3) ref+=field;
    else if (field_id==4) {
        alt+=field;
        is_symbolic=field.starts_with(SYMBOLIC_CHAR_OPEN);
        if (field.starts_with(VCF_MISSING_CHAR) || n_alts>1 || field=="<"+INS_STR+">" || field=="<"+INS_ME_STR+">" || field=="<"+INS_NOVEL_STR+">") {
            stream.ignore(numeric_limits<streamsize>::max(),LINE_END);
            return true;
        }
    }
    else if (field_id==5) {
        if (field.starts_with(VCF_MISSING_CHAR)) { qual=-1; is_high_quality=false; }
        else {
            qual=stof(field);
            is_high_quality=qual>=min_qual;
            if (high_qual_only && !is_high_quality) {
                stream.ignore(numeric_limits<streamsize>::max(),LINE_END);
                return true;
            }
        }
    }
    else if (field_id==6) {
        filter+=field;
        if (field.starts_with(VCF_MISSING_CHAR)) is_pass=false;
        else {
            is_pass=filter==PASS_STR;
            if (pass_only && !is_pass) {
                stream.ignore(numeric_limits<streamsize>::max(),LINE_END);
                return true;
            }
        }
    }
    else if (field_id==7) {
        info+=field;
        if (field.starts_with(VCF_MISSING_CHAR)) {
            sv_type=-1; sv_length=UINT_MAX;
            stream.ignore(numeric_limits<streamsize>::max(),LINE_END);
            return true;
        }
        set_sv_type(tmp_buffer);
        set_sv_length(tmp_buffer);
        if (sv_type==-1 || sv_length<min_sv_length) {
            stream.ignore(numeric_limits<streamsize>::max(),LINE_END);
            return true;
        }
    }
    else if (field_id==8) {
        format+=field;
        if (skip_samples) {
            stream.ignore(numeric_limits<streamsize>::max(),LINE_END);
            return true;
        }
    }
    else {
        n_samples++;
        if (n_samples>N_SAMPLES_IN_VCF) {
            stream.ignore(numeric_limits<streamsize>::max(),LINE_END);
            return true;
        }
        else {
            genotypes[n_samples-1].clear();
            genotypes[n_samples-1]+=field;
            ncalls_in_sample(field,is_diploid,tmp_pair);
            n_haplotypes_ref+=tmp_pair.first; n_haplotypes_alt+=tmp_pair.second;
        }
    }
    return false;
}


void VcfRecord::set_sv_type(string& tmp_buffer) {
    tmp_buffer.clear();
    const bool found = get_info_field(info,SVTYPE_STR,tmp_buffer);
    if (found) sv_type=string_to_svtype_info(tmp_buffer);
    else if (alt.starts_with('<') && alt.ends_with('>')) sv_type=string_to_svtype_alt(alt);
    else if (alt.starts_with('[') || alt.starts_with(']') || alt.ends_with('[') || alt.ends_with(']')) sv_type=string_to_svtype_alt(alt);
    else {
        const size_t ref_length = ref.length();
        const size_t alt_length = alt.length();
        if (ref_length==1 && alt_length>ref_length) sv_type=TYPE_INSERTION;
        else if (alt_length==1 && ref_length>alt_length) sv_type=TYPE_DELETION;
        else if (ref_length>1 && alt_length>1) sv_type=TYPE_REPLACEMENT;
        else sv_type=-1;
    }
}


void VcfRecord::set_sv_length(string& tmp_buffer) {
    tmp_buffer.clear();
    if (get_info_field(info,SVLEN_STR,tmp_buffer)) {
        const int32_t length = stoi(tmp_buffer);
        sv_length=length<0?-length:length;
    }
    else if (get_info_field(info,END_STR,tmp_buffer)) {
        const uint64_t end = stoi(tmp_buffer);
        sv_length=end-pos;
    }
    else if (alt.starts_with('<') && alt.ends_with('>')) sv_length=UINT_MAX;
    else if (alt.starts_with('[') || alt.starts_with(']') || alt.ends_with('[') || alt.ends_with(']')) sv_length=UINT_MAX;
    else {
        const size_t ref_length = ref.length();
        const size_t alt_length = alt.length();
        if (ref_length==1 && alt_length>ref_length) sv_length=(uint32_t)(alt_length-1);
        else if (alt_length==1 && ref_length>alt_length) sv_length=(uint32_t)(ref_length-1);
        else sv_length=UINT_MAX;
    }
}


void VcfRecord::print(ostream& stream) const {
    stream << chrom << VCF_SEPARATOR << pos << VCF_SEPARATOR << id << VCF_SEPARATOR << ref << VCF_SEPARATOR << alt << VCF_SEPARATOR;
    if (qual==-1) stream << '.';
    else stream << qual;
    stream << VCF_SEPARATOR << filter << VCF_SEPARATOR << info << VCF_SEPARATOR << format;
    for (uint32_t i=0; i<n_samples; i++) stream << VCF_SEPARATOR << genotypes[i];
}


uint8_t VcfRecord::get_gt(uint32_t sample, pair<int8_t,int8_t>& out, string& tmp_buffer) const {
    if (sample>=n_samples) return 0;
    char c;
    const size_t LENGTH = genotypes[sample].length();
    size_t i;

    out.first=-1; out.second=-1;
    if (!is_diploid) {
        if (genotypes[sample].starts_with(VCF_MISSING_CHAR)) return 1;
        tmp_buffer.clear();
        for (i=0; i<LENGTH; i++) {
            c=genotypes[sample].at(i);
            if (c==GT_SEPARATOR) {
                out.first=(int8_t)stoi(tmp_buffer);
                return 1;
            }
            else tmp_buffer+=c;
        }
        out.first=(int8_t)stoi(genotypes[sample]);
        return 1;
    }
    else {
        tmp_buffer.clear();
        for (i=0; i<LENGTH; i++) {
            c=genotypes[sample].at(i);
            if (c==PHASED_CHAR || c==UNPHASED_CHAR) {
                if (!tmp_buffer.starts_with(VCF_MISSING_CHAR)) out.first=(int8_t)stoi(tmp_buffer);
                tmp_buffer.clear();
            }
            else if (c==GT_SEPARATOR) {
                if (!tmp_buffer.starts_with(VCF_MISSING_CHAR)) out.second=(int8_t)stoi(tmp_buffer);
                return 2;
            }
            else tmp_buffer+=c;
        }
        if (!tmp_buffer.starts_with(VCF_MISSING_CHAR)) out.second=(int8_t)stoi(tmp_buffer);
        if (out.first==-1 && out.second!=-1) { out.first=out.second; out.second=-1; return 1; }
        else return 2;
    }
}


bool VcfRecord::passes_constraints() const {
    return n_alts==1 &&
           (is_high_quality || !HIGH_QUAL_ONLY) &&
           (!PASS_ONLY || is_pass || filter.starts_with(VCF_MISSING_CHAR)) &&
           sv_type!=-1 && sv_length>=MIN_SV_LENGTH &&
           (SKIP_SAMPLES || (n_samples==N_SAMPLES_IN_VCF && n_haplotypes_alt>=MIN_N_HAPLOTYPES_ALT && n_haplotypes_ref+n_haplotypes_alt>=MIN_N_HAPLOTYPES_NONMISSING));
}


bool get_info_field(const string &info, const string &key, string &out) {
    const uint8_t KEY_LENGTH = key.length()+1;
    size_t i, p, q;

    p=info.find(key+"=",0);
    if (p==string::npos) return false;
    if (key==END_STR && p>=2 && info.at(p-2)=='C' && info.at(p-1)=='I') {
        p=info.find(key+"=",p+KEY_LENGTH);
        if (p==string::npos) return false;
    }
    q=info.find(INFO_SEPARATOR,p+KEY_LENGTH);
    if (q==string::npos) q=info.length();
    out.clear();
    for (i=p+KEY_LENGTH; i<q; i++) out+=info.at(i);
    return true;
}


void ncalls_in_sample(const string& buffer, bool is_diploid, pair<uint8_t, uint8_t>& out) {
    char c;
    const size_t LENGTH = buffer.length();
    size_t i, j;

    out.first=0; out.second=0;
    if (!is_diploid) {
        if (buffer.starts_with('0')) out.first=1;
        else if (!buffer.starts_with(VCF_MISSING_CHAR)) out.second=1;
    }
    else {
        j=0;
        for (i=0; i<LENGTH; i++) {
            c=buffer.at(i);
            if (c==PHASED_CHAR || c==UNPHASED_CHAR) {
                if (buffer.starts_with('0')) out.first=1;
                else if (!buffer.starts_with(VCF_MISSING_CHAR)) out.second=1;
                j=i+1;
            }
            else if (c==GT_SEPARATOR) {
                c=buffer.at(j);
                if (c=='0') out.first++;
                else if (c!=VCF_MISSING_CHAR) out.second++;
                return;
            }
        }
        c=buffer.at(j);
        if (c=='0') out.first++;
        else if (c!=VCF_MISSING_CHAR) out.second++;
    }
}


bool is_phased(const string& buffer) {
    const size_t LENGTH = buffer.length();

    for (size_t i=0; i<LENGTH; i++) {
        switch (buffer.at(i)) {
            case PHASED_CHAR: return true;
            case UNPHASED_CHAR: return false;
            case GT_SEPARATOR: break;
        }
    }
    return true;
}


void info2map(const string& info, unordered_map<string,string>& map, string& tmp_buffer_1, string& tmp_buffer_2) {
    bool found;
    char c;
    uint16_t i;
    const uint16_t LENGTH = info.length();

    map.clear(); found=false;
    tmp_buffer_1.clear(); tmp_buffer_2.clear();
    for (i=0; i<LENGTH; i++) {
        c=info.at(i);
        if (c==INFO_ASSIGNMENT) found=true;
        else if (c==INFO_SEPARATOR) {
            map.emplace(tmp_buffer_1,tmp_buffer_2);
            tmp_buffer_1.clear(); tmp_buffer_2.clear();
            found=false;
        }
        else if (!found) tmp_buffer_1+=c;
        else tmp_buffer_2+=c;
    }
    if (found) map.emplace(tmp_buffer_1,tmp_buffer_2);
}


void for_record_in_vcf(const string& path, const function<void(VcfRecord& record)>& f, uint32_t progress_n_lines, const string& chromosome, bool is_diploid, bool high_qual_only, float min_qual, bool pass_only, uint32_t min_sv_length, bool skip_samples, uint32_t n_samples_in_vcf, float min_allele_frequency, float min_nonmissing_frequency) {
    char c;
    uint64_t n_lines;
    string tmp_buffer_1, tmp_buffer_2;
    pair<uint8_t, uint8_t> tmp_pair;
    VcfRecord record(chromosome,is_diploid,high_qual_only,min_qual,pass_only,min_sv_length,skip_samples,n_samples_in_vcf,min_allele_frequency,min_nonmissing_frequency);
    ifstream file;

    file.open(path,std::ifstream::in);
    if (!file.is_open() || !file.good()) throw runtime_error("ERROR: could not read file: " + path);
    n_lines=0;
    while (!file.eof()) {
        c=file.peek();
        if (c==VCF_COMMENT) {
            file.ignore(numeric_limits<streamsize>::max(),LINE_END);
            continue;
        }
        else if (c==EOF) break;
        record.set(file,tmp_buffer_1,tmp_buffer_2,tmp_pair);
        if (record.passes_constraints()) f(record);
        n_lines++;
        if (progress_n_lines!=0 && n_lines%progress_n_lines==0) cerr << "Scanned " << n_lines << " lines\n";
    }
    file.close();
    if (progress_n_lines!=0) cerr << "Scanned " << n_lines << " total lines\n";
}

}