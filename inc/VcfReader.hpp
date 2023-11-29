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
using std::string_view;
using std::vector;
using std::unordered_map;
using std::function;
using std::ifstream;
using std::ostream;
using std::cerr;
using std::runtime_error;

namespace sv_merge {

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


class VcfReader {
public:
    /**
     * Basic VCF constants
     */
    static const char VCF_COMMENT;
    static const char LINE_END;
    static const char VCF_SEPARATOR;
    static const char VCF_MISSING_CHAR;
    static const char SYMBOLIC_CHAR_OPEN;
    static const char SYMBOLIC_CHAR_CLOSE;
    static const string PASS_STR;
    static const char ALT_SEPARATOR;
    static const char GT_SEPARATOR;
    static const char UNPHASED_CHAR;
    static const char PHASED_CHAR;

    /**
     * Info field constants
     */
    static const char INFO_ASSIGNMENT;
    static const char INFO_SEPARATOR;
    static const string SVTYPE_STR;
    static const string SVLEN_STR;
    static const string END_STR;
    static const string CIPOS_STR;
    static const uint16_t CIPOS_STR_LENGTH;
    static const string CIEND_STR;
    static const uint16_t CIEND_STR_LENGTH;
    static const string CILEN_STR;
    static const uint16_t CILEN_STR_LENGTH;
    static const string PRECISE_STR;
    static const string IMPRECISE_STR;

    /**
     * Supported SV types
     */
    static const uint8_t TYPE_INSERTION;
    static const uint8_t TYPE_DELETION;
    static const uint8_t TYPE_INVERSION;
    static const uint8_t TYPE_DUPLICATION;
    static const uint8_t TYPE_BREAKEND;
    static const uint8_t TYPE_REPLACEMENT;

    /**
     * Supported SV types: labels used by the callers.
     */
    static const string DEL_STR;
    static const string DEL_ME_STR;
    static const string INS_STR;
    static const string INS_ME_STR;
    static const string INS_NOVEL_STR;
    static const string DUP_STR;
    static const string DUP_TANDEM_STR;
    static const string DUP_INT_STR;
    static const string INV_STR;
    static const string BND_STR;

    /**
     * Internal state of VcfReader. See `VcfRecord` for details.
     */
    uint32_t progress_n_lines;
    bool high_qual_only;
    float min_qual;
    bool pass_only;
    uint32_t min_sv_length;
    bool skip_samples;
    uint32_t n_samples_in_vcf;
    float min_allele_frequency;
    float min_nonmissing_frequency;

    /**
     * @param path a VCF file that contains only calls in `chromosome`;
     * @param callback called on every VCF record that passes the constraints; see `VcfRecord` for details;
     * @param progress_n_lines prints a progress message to STDERR after this number of lines have been read (0=silent).
     */
    VcfReader(const string& path, const function<void(VcfRecord& record)>& callback, uint32_t progress_n_lines, const string& chromosome, bool is_diploid, bool high_qual_only, float min_qual, bool pass_only, uint32_t min_sv_length, bool skip_samples, uint32_t n_samples_in_vcf, float min_allele_frequency, float min_nonmissing_frequency);
    VcfReader(const string& path, const function<void(VcfRecord& record)>& callback, const string& chromosome, bool is_diploid);

    void for_record_in_vcf();

private:
    /**
     * Internal state of VcfReader
     */
    string path;
    function<void(VcfRecord& record)> callback;
    string chromosome;
    bool is_diploid;

    /**
     * Remark: this performs a linear scan of `info`.
     *
     * @return TRUE iff `key` is found in `ìnfo`; in this case `out` contains the value of `key`.
     */
    static bool get_info_field(const string& info, const string& key, string& out);

    /**
     * Remark: the procedure handles the case where there are two copies of the chromosome but the GT field contains
     * just one value (this can happen with CNVs or when the sample is a haploid reference).
     *
     * @param buffer the content of a VCF sample field;
     * @return out output pair, containing the number of zero (`.first`) and nonzero (`.second`) haplotypes.
     */
    static void ncalls_in_sample(const string& buffer, bool is_diploid, pair<uint8_t, uint8_t>& out);

    /**
     * @param buffer the content of a VCF sample field (can contain just one haplotype).
     */
    static bool is_phased(const string& buffer);

    /**
     * Stores in `map` the set of all and only the key->value relationships in `info`.
     *
     * @param tmp_buffer_* reused temporary space.
     */
    static void info2map(const string& info, unordered_map<string,string>& map, string& tmp_buffer_1, string& tmp_buffer_2);
}

}