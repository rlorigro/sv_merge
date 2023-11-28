#pragma once

#include <climits>
#include <utility>
#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <ostream>
#include <iostream>

using std::pair;
using std::string;
using std::vector;
using std::unordered_map;
using std::function;
using std::ifstream;
using std::ostream;
using std::cerr;
using std::runtime_error;

namespace sv_merge {


/**
 * Basic VCF constants
 */
const char VCF_COMMENT = '#';
const char LINE_END = '\n';  // Change to '\r\n' for CR+LF.
const char VCF_SEPARATOR = '\t';
const char VCF_MISSING_CHAR = '.';
const char SYMBOLIC_CHAR_OPEN = '<';
const char SYMBOLIC_CHAR_CLOSE = '>';
const string PASS_STR = "PASS";
const char ALT_SEPARATOR = ',';
const char GT_SEPARATOR = ':';
const char UNPHASED_CHAR = '/';
const char PHASED_CHAR = '|';

/**
 * Info field constants
 */
const char INFO_ASSIGNMENT = '=';
const char INFO_SEPARATOR = ';';
const string SVTYPE_STR = "SVTYPE";
const string SVLEN_STR = "SVLEN";
const string END_STR = "END";
const string CIPOS_STR = "CIPOS";
const uint16_t CIPOS_STR_LENGTH = CIPOS_STR.length();
const string CIEND_STR = "CIEND";
const uint16_t CIEND_STR_LENGTH = CIEND_STR.length();
const string CILEN_STR = "CILEN";
const uint16_t CILEN_STR_LENGTH = CILEN_STR.length();
const string PRECISE_STR = "PRECISE";
const string IMPRECISE_STR = "IMPRECISE";

/**
 * Supported SV types
 */
const uint8_t TYPE_INSERTION = 1;
const uint8_t TYPE_DELETION = 2;
const uint8_t TYPE_INVERSION = 3;
const uint8_t TYPE_DUPLICATION = 4;
const uint8_t TYPE_BREAKEND = 5;
const uint8_t TYPE_REPLACEMENT = 6;

/**
 * Supported SV types: labels used by the callers.
 */
const string DEL_STR = "DEL";
const string DEL_ME_STR = "DEL:ME";
const string INS_STR = "INS";
const string INS_ME_STR = "INS:ME";
const string INS_NOVEL_STR = "INS:NOVEL";
const string DUP_STR = "DUP";
const string DUP_TANDEM_STR = "DUP:TANDEM";
const string DUP_INT_STR = "DUP:INT";
const string INV_STR = "INV";
const string BND_STR = "BND";


/**
 * Reused container, allocated only once and overwritten with each VCF line.
 *
 * Remark: we assume that every line of the VCF file belongs to the same chromosome.
 */
class VcfRecord {
public:
    /**
     * Properties that are set by the user
     */
    bool is_diploid;
    string chrom;

    /**
     * VCF fields loaded from a line
     */
    float qual;  // -1 = missing
    uint64_t pos;
    string id, ref, alt, filter, info, format;
    vector<string> genotypes;

    /**
     * Properties that are already known after the first scan of a VCF line
     */
    bool is_high_quality, is_pass, is_symbolic;
    int8_t sv_type;  // -1 = unsupported type
    uint32_t sv_length;  // UINT_MAX = the length of the SV could not be inferred
    uint32_t n_alts;  // >1 iff the site is multiallelic
    uint32_t n_samples;  // Actual number of samples in the record (if >N_SAMPLES_IN_VCF, it is fixed to N_SAMPLES_IN_VCF+1).
    uint32_t n_haplotypes_ref, n_haplotypes_alt;

    /**
     * @param chromosome the chromosome to which all lines in the VCF belong;
     * @param is_diploid TRUE iff there are two copies of `chromosome`;
     * @param n_samples_in_vcf expected number of samples in every line;
     * @param pass_only calls with FILTER different from PASS or `.` are discarded by the user;
     * @param high_qual_only calls with `QUAL<min_qual` are discarded by the user;
     * @param skip_samples sample GT fields are discarded by the user;
     * @param min_sv_length calls shorter than this are discarded by the user;
     * @param min_allele_frequency calls with too few haplotypes containing the ALT allele are discarded by the user;
     * the procedure assumes that every sample has the same ploidy, so it might not work for CNVs or if a sample is a
     * haploid reference;
     * @param min_nonmissing_frequency calls with too few non-missing haplotypes are discarded by the user;
     * the procedure assumes that every sample has the same ploidy, so it might not work for CNVs or if a sample is a
     * haploid reference.
     */
    VcfRecord(const string& chromosome, bool is_diploid, bool high_qual_only, float min_qual, bool pass_only, uint32_t min_sv_length, bool skip_samples, uint32_t n_samples_in_vcf, float min_allele_frequency, float min_nonmissing_frequency);

    /**
     * Reads `stream` until EOL/EOF and loads some or all of the data in the current line.
     *
     * - If the call is a symbolic INS, no field after ALT is loaded.
     * - If there are zero or more than one ALT, no field after ALT is loaded.
     * - If `HIGH_QUAL_ONLY` is true and the call has low quality, no field after QUAL is loaded.
     * - If `PASS_ONLY` is true and the call is not PASS, no field after FILTER is loaded.
     * - If the call is not of a supported type, no field after INFO is loaded.
     * - If the call is shorter than `MIN_SV_LENGTH`, no field after INFO is loaded.
     * - If there are more samples than `N_SAMPLES_IN_VCF`, no sample after the maximum is loaded.
     *
     * @param tmp_* reused temporary space.
     */
    void set(ifstream& stream, string& tmp_buffer_1, string& tmp_buffer_2, pair<uint8_t, uint8_t>& tmp_pair);

    /**
     * @return TRUE iff the record passes all the constraints set at construction time.
     */
    bool passes_constraints() const;

    /**
     * Appends to `stream` a string representation of the object. No terminator is added at the end.
     */
    void print(ostream& stream) const;

    /**
     * Remark: the procedure handles the case where there are two copies of the chromosome but the GT field contains
     * just one value (this can happen with CNVs or when the sample is a haploid reference).
     *
     * @param sample sample ID (the first sample has ID zero);
     * @param out output pair that stores the only GT haplotype ID in `first`, or the two haplotype IDs in `first` and
     * `second`, depending on the ploidy of `chrom`; -1 = missing haplotype;
     * @param tmp_buffer reused temporary space;
     * @return number of GT haplotypes in `sample` (zero iff `sample` is invalid).
     */
    uint8_t get_gt(uint32_t sample, pair<int8_t,int8_t>& out, string& tmp_buffer) const;

private:
    /**
     * Constraints specified by the user at construction time.
     *
     * Remark: PASS_ONLY means FILTER=PASS or FILTER='.'
     */
    bool PASS_ONLY, HIGH_QUAL_ONLY, SKIP_SAMPLES;
    uint32_t MIN_SV_LENGTH, N_SAMPLES_IN_VCF;
    float MIN_QUAL;
    double MIN_N_HAPLOTYPES_ALT, MIN_N_HAPLOTYPES_NONMISSING;

    /**
     * Sets `sv_type` using `ref`, `alt` and `info`, which are assumed to be already set.
     *
     * @param tmp_buffer reused temporary space.
     */
    void set_sv_type(string& tmp_buffer);

    /**
     * Sets `sv_length` using the following VCF fields, in order:
     * - INFO.SVLEN if it exists;
     * - INFO.END if it exists;
     * - non-symbolic ALT if it exists.
     *
     * Remark: `ref`, `alt` and `info` are assumed to be already set.
     *
     * @param tmp_buffer reused temporary space.
     */
    void set_sv_length(string& tmp_buffer);

    /**
     * Main logic of `set()`.
     *
     * @param tmp_* reused temporary space;
     * @return TRUE iff some VCF fields were skipped.
     */
    bool set_field(const string& field, uint32_t field_id, bool high_qual_only, float min_qual, bool pass_only, uint32_t min_sv_length, bool skip_samples, ifstream& stream, string& tmp_buffer, pair<uint8_t, uint8_t>& tmp_pair);
};

/**
 * Calls `f` on every VCF record that passes the constraints. See `VcfRecord` for details.
 *
 * @param path a VCF file that contains only calls in `chromosome`;
 * @param progress_n_lines prints a progress message to STDERR after this number of lines have been read (0=silent).
 */
void for_record_in_vcf(const string& path, const function<void(VcfRecord& record)>& f, uint32_t progress_n_lines, const string& chromosome, bool is_diploid, bool high_qual_only, float min_qual, bool pass_only, uint32_t min_sv_length, bool skip_samples, uint32_t n_samples_in_vcf, float min_allele_frequency, float min_nonmissing_frequency);

/**
 * Remark: this performs a linear scan of `info`.
 *
 * @return TRUE iff `key` is found in `ìnfo`; in this case `out` contains the value of `key`.
 */
bool get_info_field(const string& info, const string& key, string& out);

/**
 * Remark: the procedure handles the case where there are two copies of the chromosome but the GT field contains
 * just one value (this can happen with CNVs or when the sample is a haploid reference).
 *
 * @param buffer the content of a VCF sample field;
 * @return out output pair, containing the number of zero (`.first`) and nonzero (`.second`) haplotypes.
 */
void ncalls_in_sample(const string& buffer, bool is_diploid, pair<uint8_t, uint8_t>& out);

/**
 * @param buffer the content of a VCF sample field (can contain just one haplotype).
 */
bool is_phased(const string& buffer);

/**
 * Stores in `map` the set of all and only the key->value relationships in `info`.
 *
 * @param tmp_buffer_* reused temporary space.
 */
void info2map(const string& info, unordered_map<string,string>& map, string& tmp_buffer_1, string& tmp_buffer_2);

}