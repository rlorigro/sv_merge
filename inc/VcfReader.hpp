#pragma once

#include "Filesystem.hpp"

using ghc::filesystem::path;

#include <climits>
#include <numeric>
#include <cmath>
#include <utility>
#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <functional>
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
using std::ceil;
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
     * Max size of a file stream
     */
    static const uint64_t STREAMSIZE_MAX;

    /**
     * Constraints specified by the user at construction time.
     *
     * Remark: `pass_only` means FILTER=PASS or FILTER='.'
     */
    bool pass_only, high_qual_only;
    uint32_t min_sv_length, n_samples_to_load;
    float min_qual, min_af, min_nmf;

    /**
     * Values derived from the user constraints above
     */
    float min_n_haplotypes_alt_baseline, min_n_haplotypes_nonmissing_baseline;
    uint32_t min_n_haplotypes_alt, min_n_haplotypes_nonmissing;

    /**
     * VCF fields loaded from a line
     */
    bool is_autosomal;
    float qual;  // -1 = missing
    uint64_t pos;
    string chrom, id, ref, alt, filter, info, format;
    vector<string> genotypes;

    /**
     * Properties that are already known after the first scan of a VCF line
     */
    bool is_high_quality, is_pass, is_symbolic;
    int8_t sv_type;  // -1 = unsupported type
    uint32_t sv_length;  // MAX = the length of the SV could not be inferred
    uint32_t n_alts;  // >1 iff the site is multiallelic
    uint32_t n_samples;  // Actual number of samples in the record (if >n_samples_to_load, it is fixed to n_samples_to_load+1).
    uint32_t n_haplotypes_ref, n_haplotypes_alt;

    /**
     * @param n_samples_to_load number of samples in a line that are interesting for the user;
     * @param pass_only calls with FILTER different from PASS or `.` are discarded by the user;
     * @param high_qual_only calls with `QUAL<min_qual` are discarded by the user;
     * @param min_sv_length calls shorter than this are discarded by the user;
     * @param min_af calls with too few haplotypes containing the ALT allele are discarded by the user;
     * @param min_nmf calls with too few non-missing haplotypes are discarded by the user.
     */
    VcfRecord(bool high_qual_only, float min_qual, bool pass_only, uint32_t min_sv_length, uint32_t n_samples_to_load, float min_af, float min_nmf);

    /**
     * @return a new object that contains a copy of every variable of the current object that was loaded from a VCF
     * record; the state of every other variable is undefined.
     */
    VcfRecord clone() const;

    /**
     * Reads `stream` until EOL/EOF and loads some or all of the data in the current line.
     *
     * - If the call is a symbolic INS, no field after ALT is loaded.
     * - If there are zero or more than one ALT, no field after ALT is loaded.
     * - If `high_qual_only` is true and the call has low quality, no field after QUAL is loaded.
     * - If `pass_only` is true and the call is not PASS, no field after FILTER is loaded.
     * - If the call is not of a supported type, no field after INFO is loaded.
     * - If the call is shorter than `min_sv_length`, no field after INFO is loaded.
     * - If there are more samples than `n_samples_to_load`, no sample after the maximum is loaded.
     */
    void set(ifstream& stream);

    /**
     * @return TRUE iff the record passes all the user constraints set at construction time.
     */
    bool passes_constraints() const;

    /**
     * Appends to `stream` a string representation of the object. No terminator is added at the end.
     */
    void print(ostream& stream) const;

    /**
     * Remark: it can happen that there are two copies of a chromosome but the GT field contains just one value (e.g.
     * with CNVs).
     *
     * @param sample sample ID (the first sample has ID zero);
     * @param out output pair that stores the only GT haplotype ID in `first`, or the two haplotype IDs in `first` and
     * `second`, depending on the ploidy of `chrom` in `sample`; -1 = missing haplotype;
     * @param tmp_buffer reused temporary space;
     * @return number of GT haplotypes in `sample` (zero iff `sample` is invalid).
     */
    uint8_t get_gt(uint32_t sample, pair<int8_t,int8_t>& out);

    /**
     * Remark: this performs a linear scan of `info`.
     *
     * @return the first occurrence of `key` in `ìnfo[from..]`, if present (in this case `out` contains the value of
     * `key`); `string::npos` if `key` does not occur in `ìnfo`.
     */
    size_t get_info_field(const string& key, const size_t from, string& out) const;

    /**
     * Remark: it can happen that there are two copies of a chromosome but the GT field contains just one value (e.g.
     * with CNVs).
     *
     * @param buffer the content of a VCF sample field;
     * @return out output pair, containing the number of zero (`.first`) and nonzero (`.second`) haplotypes.
     */
    void ncalls_in_sample(const string& buffer, pair<uint8_t, uint8_t>& out) const;

    /**
     * Remark: a haploid sample is assumed to be phased.
     *
     * @param buffer the content of a VCF sample field (can contain just one haplotype).
     */
    static bool is_phased(const string& buffer);

    /**
     * Stores in `map` the set of all and only the key->value relationships in `info`.
     */
    void info2map(unordered_map<string,string>& map);

    /**
     * Stores in `out` the zero-based IDs of all the elements of `genotypes` that contain a nonzero.
     */
    void get_samples_with_alt(vector<uint32_t>& out);

    /**
     * Same as above, but returns just the string ID of each sample.
     */
    void get_samples_with_alt(vector<string>& out);

    bool is_alt_symbolic() const;

    /**
     * @return 0=single breakend without inserted sequence; 1=with inserted sequence; 2=not a single breakend.
     */
    uint8_t is_breakend_single() const;

    /**
     * Virtual telomeric breakends are artificial records that carry no information.
     */
    bool is_breakend_virtual(const unordered_map<string,string>& chromosomes);

    /**
     * Remark: `out` is set to an empty string if the chromosome could not be determined.
     */
    void get_breakend_chromosome(string& out) const;

    /**
     * @return UINT32_MAX if the position could not be determined; otherwise, the original, one-based value.
     */
    uint32_t get_breakend_pos();

    /**
     * Remark: the procedure works also for single breakends.
     *
     * @return
     * 0 if the orientation could not be determined;
     * 1 if the CIS side of the breakend extends to the left of `pos`;
     * 2 if the CIS side of the breakend extends the the right of `pos`.
     */
    uint8_t get_breakend_orientation_cis() const;

    /**
    * @return
    * 0 if the orientation could not be determined;
    * 1 if the TRANS side of the breakend extends to the left of the position in `alt`;
    * 2 if the TRANS side of the breakend extends the the right of the position in `alt`.
    */
    uint8_t get_breakend_orientation_trans() const;

    /**
     * Stores in `out` all bases inserted between the breakend position and its mate. I.e. if REF=`T` and
     * ALT=`]chr13:123456]AGTNNNNNCAT`, the procedure returns `AGTNNNNNCA`.
     *
     * Remark: the procedure works also for single breakends.
     */
    void get_breakend_inserted_sequence(string& out) const;

    /**
     * Checks the IMPRECISE tag and the confidence intervals fields.
     */
    bool is_precise();

    /**
	 * @param out `first`: quantity to be added to `pos` to get the first position of the confidence interval
     * (typically <=0); `second`: quantity to be added to `pos` to get the last position of the confidence interval
     * (typically >=0).
	 */
    void get_confidence_interval_pos(pair<float, float>& out);
    void get_confidence_interval_end(pair<float, float>& out);
    void get_confidence_interval_length(pair<float, float>& out);

    /**
     * Stores in `out` the zero-based coordinates of the SV in the reference.
     * - If the SV affects all and only the zero-based positions in a closed interval `[x..y]`, `out=(x,y+1)`.
     * - If the SV is an INS between zero-based positions `x` and `x+1`, `out=(x+1,x+1)`. Note that `x=-1` is allowed
     *   (telomeric insertion).
     * - If the SV is a BND, either between zero-based positions `x` and `x+1`, or between zero-based positions `x-1`
     *   and `x`, `out=(x,x)`. The spec allows `x=-1` (virtual telomeric breakend), but the function returns UINT32_MAX
     *   instead, since a virtual breakend carries no information.
     * - If a value cannot be determined, it is set to UINT32_MAX.
     *
     * @param use_confidence_intervals enlarges the coordinates above using confidence intervals on `pos` and
     * `sv_length`, if available.
     */
    void get_reference_coordinates(bool use_confidence_intervals, pair<uint32_t, uint32_t>& out);

private:
    /**
     * Reused temporary space
     */
    string tmp_buffer_1, tmp_buffer_2;
    pair<uint8_t, uint8_t> tmp_pair;
    pair<float, float> tmp_pair_2;

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
     * Remark: the length of a replacement record is arbitrarily set to the length of the substring of the reference
     * that is to be replaced.
     * Remark: `ref`, `alt` and `info` are assumed to be already set.
     *
     * @param tmp_buffer reused temporary space.
     */
    void set_sv_length(string& tmp_buffer);

    /**
     * Main logic of `set()`.
     *
     * @return TRUE iff some VCF fields were skipped.
     */
    bool set_field(const string& field, uint32_t field_id, bool high_qual_only, float min_qual, bool pass_only, uint32_t min_sv_length, float min_af, float min_nmf, ifstream& stream, string& tmp_buffer, pair<uint8_t, uint8_t>& tmp_pair);

    /**
     * Core logic of confidence intervals extraction
     *
     * @param which 0=pos, 1=sv_length, 2=end.
     */
    void get_confidence_interval(uint8_t which, pair<float, float>& out);
};


class VcfReader {
public:
    /**
     * Basic VCF constants
     */
    static const char VCF_COMMENT;
    static const string VCF_HEADER_PREFIX;
    static const uint8_t VCF_HEADER_PREFIX_LENGTH;
    static const uint8_t N_NONSAMPLE_FIELDS_VCF;  // Including FORMAT
    static const char LINE_END;
    static const char VCF_SEPARATOR;
    static const char VCF_MISSING_CHAR;
    static const char SYMBOLIC_CHAR_OPEN;
    static const char SYMBOLIC_CHAR_CLOSE;
    static const char BREAKEND_CHAR_OPEN;
    static const char BREAKEND_CHAR_CLOSE;
    static const string PASS_STR;
    static const char ALT_SEPARATOR;
    static const char GT_SEPARATOR;
    static const char BREAKEND_SEPARATOR;
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
    static const char CI_SEPARATOR;
    static const string STDEV_POS_STR;
    static const uint16_t STDEV_POS_STR_LENGTH;
    static const string STDEV_LEN_STR;
    static const uint16_t STDEV_LEN_STR_LENGTH;
    static const string STDEV_END_STR;
    static const uint16_t STDEV_END_STR_LENGTH;
    static const string PRECISE_STR;
    static const string IMPRECISE_STR;
    static const string MATEID_STR;

    /**
     * Supported SV types
     */
    static const uint8_t TYPE_INSERTION;
    static const uint8_t TYPE_DELETION;
    static const uint8_t TYPE_INVERSION;
    static const uint8_t TYPE_DUPLICATION;
    static const uint8_t TYPE_BREAKEND;
    static const uint8_t TYPE_REPLACEMENT;
    static const uint8_t TYPE_CNV;

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
    static const string CNV_STR;

    /**
     * Configuration parameters. See `VcfRecord` for details.
     * This is public to allow users to set only specific parameters.
     */
    uint32_t progress_n_lines;
    bool high_qual_only;
    float min_qual;
    bool pass_only;
    uint32_t min_sv_length;
    uint32_t n_samples_to_load;  // A prefix of the list of all samples. 0=do not load any sample. MAX=load all samples.
    float min_allele_frequency;
    float min_nonmissing_frequency;

    /**
     * All the samples in the VCF according to the header
     */
    vector<string> sample_ids;
    uint32_t n_samples_in_vcf;

    /**
     * @param path a VCF file that contains only calls in `chromosome`;
     * @param callback called on every VCF record that passes the constraints; see `VcfRecord` for details;
     * @param progress_n_lines prints a progress message to STDERR after this number of lines have been read (0=silent).
     */
    VcfReader(const path& vcf_path, uint32_t progress_n_lines, bool high_qual_only, float min_qual, bool pass_only, uint32_t min_sv_length, uint32_t n_samples_to_load, float min_allele_frequency, float min_nonmissing_frequency);
    VcfReader(const path& vcf_path);

    void for_record_in_vcf(const function<void(VcfRecord& record)>& callback);

private:
    /**
     * Internal state of VcfReader
     */
    string vcf_path;
    function<void(VcfRecord& record)> callback;
};

}