#include "VcfReader.hpp"
#include "VariantGraph.hpp"
#include "misc.hpp"

using sv_merge::interval_t;
using sv_merge::run_command;
using sv_merge::VcfReader;
using sv_merge::VariantGraph;

#include <iostream>
#include <algorithm>
#include <random>

using std::ofstream;


void print_truth_vcf_header(ofstream& out) {
    out << "##fileformat=VCFv4.2\n";
    out << "##contig=<ID=chr1,length=131>\n";
    out << "##contig=<ID=chr2,length=285>\n";
    out << "##contig=<ID=chr3,length=67>\n";
    out << R"(##ALT=<ID=INS,Description="Insertion">)" << '\n';
    out << R"(##ALT=<ID=DEL,Description="Deletion">)" << '\n';
    out << R"(##ALT=<ID=DUP,Description="Duplication">)" << '\n';
    out << R"(##ALT=<ID=INV,Description="Inversion">)" << '\n';
    out << R"(##ALT=<ID=BND,Description="Breakend">)" << '\n';
    out << R"(##FILTER=<ID=PASS,Description="All filters passed">)" << '\n';
    out << R"(##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">)" << '\n';
    out << R"(##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variation">)" << '\n';
    out << R"(##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variation">)" << '\n';
    out << R"(##INFO=<ID=CHR2,Number=1,Type=String,Description="Mate chromosome for BND SVs">)" << '\n';
    out << R"(##INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variation">)" << '\n';
    out << R"(##INFO=<ID=MATEID,Number=A,Type=String,Description="ID of mate breakends">)" << '\n';
    out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
}

/**
 * @param mode 0=prints the entire input VCF; 1=only records that are loaded by `VcfReader`; 2=only records that are
 * used to build the graph.
 */
void print_truth_vcf(uint8_t mode, bool print_mateid_tag, ofstream& out) {
    const string QUAL = "60";  // Arbitrary
    const string FILTER = "PASS";
    const string FORMAT = "GT";
    const string GT = "0/1";  // Arbitrary
    const string INFIX = QUAL+"\t"+FILTER;
    const string SUFFIX = FORMAT+"\t"+GT+"\n";

    out << "chr1\t1\tsniffles_dup1\tA\t<DUP>\t" << INFIX << "\tSVTYPE=DUP;SVLEN=26;\t" << SUFFIX;  // Identical to bnd19
    out << "chr1\t2\tsniffles_bnd19\tA\t]chr1:27]A\t" << INFIX << "\tSVTYPE=BND" << (print_mateid_tag?";MATEID=bnd20\t":"\t") << SUFFIX;
    out << "chr1\t2\tsniffles_del1\tAAAAAAA\tA\t" << INFIX << "\tSVTYPE=DEL;SVLEN=6;\t" << SUFFIX;
    out << "chr1\t3\tsniffles_inv1\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=15;\t" << SUFFIX;
    out << "chr1\t3\tsniffles_bnd1\tA\t]chr2:49]A\t" << INFIX << "\tSVTYPE=BND;\t" << SUFFIX;
    if (mode<2) out << "chr1\t4\tsniffles_bnd17\tA\tA.\t" << INFIX << "\tSVTYPE=BND;\t" << SUFFIX;  // Single, skipped.
    out << "chr1\t5\tpbsv_del2\tAAAAAAAA\tA\t" << INFIX << "\tSVTYPE=DEL;SVLEN=7;\t" << SUFFIX;
    out << "chr1\t5\tpbsv_ins1\tA\tATTTTTTTTTTTTTTTTTTTT\t" << INFIX << "\tSVTYPE=INS;SVLEN=20;\t" << SUFFIX;
    if (mode<2) out << "chr1\t5\tpbsv_bnd18\tA\t.A\t" << INFIX << "\tSVTYPE=BND;\t" << SUFFIX;  // Single, skipped.
    out << "chr1\t6\tpbsv_inv2\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=15;\t" << SUFFIX;
    out << "chr1\t7\tpbsv_del3\tAAAAAAAAAA\tA\t" << INFIX << "\tSVTYPE=DEL;SVLEN=9;\t" << SUFFIX;
    out << "chr1\t7\tpbsv_bnd2\tA\tA]chr2:54]\t" << INFIX << "\tSVTYPE=BND;\t" << SUFFIX;
    out << "chr1\t8\tpbsv_rep1\tAAAAAAAAAAAAAAA\tATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\t" << INFIX << "\t.\t" << SUFFIX;  // Replacement
    out << "chr1\t8\tpav_bnd3\tA\t[chr2:58[A\t" << INFIX << "\tSVTYPE=BND;\t" << SUFFIX;
    out << "chr1\t9\tpav_ins2\tA\tATTTTTTTTTTTTTTTTTTTTTTTTT\t" << INFIX << "\tSVTYPE=INS;SVLEN=25;\t" << SUFFIX;
    out << "chr1\t10\tpav_inv3\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=14;\t" << SUFFIX;
    out << "chr1\t11\tpav_ins3\tA\tATTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\t" << INFIX << "\tSVTYPE=INS;SVLEN=30;\t" << SUFFIX;
    out << "chr1\t13\tpav_bnd9\tA\tA[chr2:167[\t" << INFIX << "\tSVTYPE=BND;\t" << SUFFIX;
    if (mode<2) out << "chr1\t15\tpav_bnd15\tA\tA.\t" << INFIX << "\tSVTYPE=BND;\t" << SUFFIX;  // Single, skipped.
    if (mode<2) out << "chr1\t16\tpav_bnd16\tA\t.A\t" << INFIX << "\tSVTYPE=BND;\t" << SUFFIX;  // Single, skipped.
    out << "chr1\t18\tsvim_bnd10\tA\t.TTTTTA\t" << INFIX << "\tSVTYPE=BND;\t" << SUFFIX;  // Single with INS
    out << "chr1\t19\tsvim_bnd4\tA\tA[chr2:61[\t" << INFIX << "\tSVTYPE=BND;\t" << SUFFIX;
    out << "chr1\t20\tsvim_inv4\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=8;\t" << SUFFIX;
    out << "chr1\t20\tsvim_bnd5\tA\t[chr2:172[TTTTTA\t" << INFIX << "\tSVTYPE=BND;\t" << SUFFIX;
    if (mode<2) out << "chr1\t23\tsvim_bnd14\tA\tA[chr3:68[\t" << INFIX << "\tSVTYPE=BND;\t" << SUFFIX;  // Virtual telomeric, skipped.
    out << "chr1\t24\tsvim_bnd6\tA\t]chr3:67]TTTTTTTTTTA\t" << INFIX << "\tSVTYPE=BND;\t" << SUFFIX;
    out << "chr1\t25\tsvim_bnd11\tA\tATTTTTTTTTTTTTTT.\t" << INFIX << "\tSVTYPE=BND;\t" << SUFFIX;  // Single with INS
    out << "chr1\t26\tcutesv_bnd12\tA\t.TTTTTTTTTTA\t" << INFIX << "\tSVTYPE=BND;\t" << SUFFIX;  // Single with INS
    out << "chr1\t27\tcutesv_bnd20\tA\tA[chr1:2[\t" << INFIX << "\tSVTYPE=BND" << (print_mateid_tag?";MATEID=bnd19\t":"\t") << SUFFIX;  // Mate, skipped.
    out << "chr1\t129\tcutesv_bnd7\tA\tATTTTTTTTTTTTTTT[chr3:60[\t" << INFIX << "\tSVTYPE=BND;\t" << SUFFIX;
    if (mode==0) out << "chr1\t130\tcutesv_ins4\tA\t<INS>\t" << INFIX << "\tSVTYPE=INS;SVLEN=30;\t" << SUFFIX;  // Symbolic INS, skipped.
    out << "chr1\t131\tcutesv_bnd8\tA\tATTTTTTTTTTTTTTTTTTTT]chr2:285]\t" << INFIX << "\tSVTYPE=BND;\t" << SUFFIX;
    if (mode<2) out << "chr1\t132\tcutesv_bnd13\tN\t[chr2:286[.\t" << INFIX << "\tSVTYPE=BND;\t" << SUFFIX;  // Virtual telomeric, skipped.

    // Additional events at the ends of chr1
    if (mode<2) {
        out << "chr1\t0\tsvision_del0\tNAAAAA\tN\t" << INFIX << '\t' << "SVTYPE=DEL\t" << SUFFIX;
        out << "chr1\t129\tsvision_del4\tAAA\tA\t" << INFIX << '\t' << "SVTYPE=DEL\t" << SUFFIX;
    }
    out << "chr1\t0\tsvision_ins0\tN\tNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\t" << INFIX << '\t' << "SVTYPE=INS\t" << SUFFIX;
    out << "chr1\t131\tsvision_ins5\tA\tATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\t" << INFIX << '\t' << "SVTYPE=INS\t" << SUFFIX;
    out << "chr1\t0\tsvision_dup2\tA\t<DUP>\t" << INFIX << '\t' << "SVTYPE=DUP;SVLEN=12\t" << SUFFIX;
    out << "chr1\t129\tsvision_dup3\tA\t<DUP>\t" << INFIX << '\t' << "SVTYPE=DUP;SVLEN=2\t" << SUFFIX;
    out << "chr1\t0\tsvision_inv5\tA\t<INV>\t" << INFIX << '\t' << "SVTYPE=INV;SVLEN=16\t" << SUFFIX;
    out << "chr1\t129\tsvision_inv6\tA\t<INV>\t" << INFIX << '\t' << "SVTYPE=INV;SVLEN=2\t" << SUFFIX;
    out << "chr1\t0\tsvision_rep2\tNAAAAAAAAAAAAAAAAA\tNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\t" << INFIX << "\t.\t" << SUFFIX;
    out << "chr1\t129\tsvision_rep3\tAAA\tATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\t" << INFIX << "\t.\t" << SUFFIX;
    out << "chr1\t1\tsvision_bnd21\tA\t]chr2:285]A\t" << INFIX << '\t' << "SVTYPE=BND\t" << SUFFIX;
}

/**
 * Excluding symbolic INS
 */
uint16_t get_n_vcf_records() { return 43; }

vector<string> get_caller_ids() { return vector<string> {"sniffles","pbsv","pav","svim","cutesv","svision"}; }

/**
 * Assumes a tandem flank equal to 10.
 *
 * @param print_paths_* prints paths that cover all SVs of a given type (possibly multiple paths per event).
 */
void print_truth_gfa(bool print_paths_del, bool print_paths_inv, bool print_paths_ins, bool print_paths_rep, bool print_paths_dup, bool print_paths_bnd, ofstream& out) {
    // chr1 ref
    out << "S\t1_1\tA\n";
    out << "S\t1_2\tA\n";
    out << "S\t1_3\tA\n";
    out << "S\t1_4\tAA\n";
    out << "S\t1_6\tA\n";
    out << "S\t1_7\tA\n";
    out << "S\t1_8\tA\n";
    out << "S\t1_9\tA\n";
    out << "S\t1_10\tA\n";
    out << "S\t1_11\tA\n";
    out << "S\t1_12\tA\n";
    out << "S\t1_13\tA\n";
    out << "S\t1_14\tAAA\n";
    out << "S\t1_17\tA\n";
    out << "S\t1_18\tA\n";
    out << "S\t1_19\tA\n";
    out << "S\t1_20\tA\n";
    out << "S\t1_21\tA\n";
    out << "S\t1_22\tA\n";
    out << "S\t1_23\tA\n";
    out << "S\t1_24\tA\n";
    out << "S\t1_25\tA\n";
    out << "S\t1_26\tAA\n";
    out << "S\t1_28\tA\n";
    out << "S\t1_29\tAAAAAAAAAAAAAAAAAAAAAAA\n";  // With tandem flank
    out << "S\t1_129\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";  // With tandem flank
    out << "S\t1_130\tAA\n";
    out << "L\t1_1\t+\t1_2\t+\t*\n";
    out << "L\t1_2\t+\t1_3\t+\t*\n";
    out << "L\t1_3\t+\t1_4\t+\t*\n";
    out << "L\t1_4\t+\t1_6\t+\t*\n";
    out << "L\t1_6\t+\t1_7\t+\t*\n";
    out << "L\t1_7\t+\t1_8\t+\t*\n";
    out << "L\t1_8\t+\t1_9\t+\t*\n";
    out << "L\t1_9\t+\t1_10\t+\t*\n";
    out << "L\t1_10\t+\t1_11\t+\t*\n";
    out << "L\t1_11\t+\t1_12\t+\t*\n";
    out << "L\t1_12\t+\t1_13\t+\t*\n";
    out << "L\t1_13\t+\t1_14\t+\t*\n";
    out << "L\t1_14\t+\t1_17\t+\t*\n";
    out << "L\t1_17\t+\t1_18\t+\t*\n";
    out << "L\t1_18\t+\t1_19\t+\t*\n";
    out << "L\t1_19\t+\t1_20\t+\t*\n";
    out << "L\t1_20\t+\t1_21\t+\t*\n";
    out << "L\t1_21\t+\t1_22\t+\t*\n";
    out << "L\t1_22\t+\t1_23\t+\t*\n";
    out << "L\t1_23\t+\t1_24\t+\t*\n";
    out << "L\t1_24\t+\t1_25\t+\t*\n";
    out << "L\t1_25\t+\t1_26\t+\t*\n";
    out << "L\t1_26\t+\t1_28\t+\t*\n";
    out << "L\t1_28\t+\t1_29\t+\t*\n";
    out << "L\t1_129\t+\t1_130\t+\t*\n";

    // chr2 ref
    out << "S\t2_40\tCCCCCCCCCC\n";
    out << "S\t2_50\tCCCCC\n";
    out << "S\t2_55\tCCC\n";
    out << "S\t2_58\tCCC\n";
    out << "S\t2_61\tCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n";  // With tandem flank
    out << "S\t2_166\tCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n";  // With tandem flank
    out << "S\t2_167\tCCCCC\n";
    out << "S\t2_172\tCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n";  // With tandem flank
    out << "S\t2_285\tCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n";  // With tandem flank
    out << "L\t2_40\t+\t2_50\t+\t*\n";
    out << "L\t2_50\t+\t2_55\t+\t*\n";
    out << "L\t2_55\t+\t2_58\t+\t*\n";
    out << "L\t2_58\t+\t2_61\t+\t*\n";
    out << "L\t2_166\t+\t2_167\t+\t*\n";
    out << "L\t2_167\t+\t2_172\t+\t*\n";

    // chr3 ref
    out << "S\t3_60\tGGGGGGGGGGGGGGGGGGGG\n";  // With tandem flank
    out << "S\t3_67\tGGGGGGGG\n";
    out << "L\t3_60\t+\t3_67\t+\t*\n";

    // chr1 non-ref
    out << "L\t1_2\t+\t1_9\t+\t*\n";  // del1
    out << "L\t1_4\t+\t1_13\t+\t*\n";  // del2
    out << "L\t1_7\t+\t1_17\t+\t*\n";  // del3

    out << "L\t1_3\t+\t1_18\t-\t*\n";  // inv1
    out << "L\t1_4\t-\t1_19\t+\t*\n";  // inv1
    out << "L\t1_6\t+\t1_21\t-\t*\n";  // inv2
    out << "L\t1_7\t-\t1_22\t+\t*\n";  // inv2
    out << "L\t1_10\t+\t1_24\t-\t*\n";  // inv3
    out << "L\t1_11\t-\t1_25\t+\t*\n";  // inv3
    out << "L\t1_20\t+\t1_28\t-\t*\n";  // inv4
    out << "L\t1_21\t-\t1_29\t+\t*\n";  // inv4

    out << "S\tins1\tTTTTTTTTTTTTTTTTTTTT\n";  // ins1
    out << "L\t1_4\t+\tins1\t+\t*\n";  // ins1
    out << "L\tins1\t+\t1_6\t+\t*\n";  // ins1
    out << "S\tins2\tTTTTTTTTTTTTTTTTTTTTTTTTT\n";  // ins2
    out << "L\t1_9\t+\tins2\t+\t*\n";  // ins2
    out << "L\tins2\t+\t1_10\t+\t*\n";  // ins2
    out << "S\tins3\tTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n";  // ins3
    out << "L\t1_11\t+\tins3\t+\t*\n";  // ins3
    out << "L\tins3\t+\t1_12\t+\t*\n";  // ins3

    out << "S\trep1\tTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n";  // rep1
    out << "L\t1_8\t+\trep1\t+\t*\n";  // rep1
    out << "L\trep1\t+\t1_23\t+\t*\n";  // rep1

    out << "L\t1_26\t+\t1_2\t+\t*\n";  // dup1 = bnd19
    out << "S\tbnd10\tTTTTT\n";  // bnd10
    out << "L\tbnd10\t+\t1_18\t+\t*\n";  // bnd10
    out << "S\tbnd11\tTTTTTTTTTTTTTTT\n";  // bnd11
    out << "L\t1_25\t+\tbnd11\t+\t*\n";  // bnd11
    out << "S\tbnd12\tTTTTTTTTTT\n";  // bnd12
    out << "L\tbnd12\t+\t1_26\t+\t*\n";  // bnd12

    out << "L\t2_40\t+\t1_3\t+\t*\n";  // bnd1
    out << "L\t1_7\t+\t2_50\t-\t*\n";  // bnd2
    out << "L\t1_8\t-\t2_58\t+\t*\n";  // bnd3
    out << "L\t1_19\t+\t2_61\t+\t*\n";  // bnd4
    out << "S\tbnd5\tTTTTT\n";  // bnd5
    out << "L\tbnd5\t+\t1_20\t+\t*\n";  // bnd5
    out << "L\t2_172\t-\tbnd5\t+\t*\n";  // bnd5
    out << "S\tbnd6\tTTTTTTTTTT\n";  // bnd6
    out << "L\t3_67\t+\tbnd6\t+\t*\n";  // bnd6
    out << "L\tbnd6\t+\t1_24\t+\t*\n";  // bnd6
    out << "S\tbnd7\tTTTTTTTTTTTTTTT\n";  // bnd7
    out << "L\t1_129\t+\tbnd7\t+\t*\n";  // bnd7
    out << "L\tbnd7\t+\t3_67\t+\t*\n";  // bnd7
    out << "S\tbnd8\tTTTTTTTTTTTTTTTTTTTT\n";  // bnd8
    out << "L\t1_130\t+\tbnd8\t+\t*\n";  // bnd8
    out << "L\tbnd8\t+\t2_285\t-\t*\n";  // bnd8
    out << "L\t1_13\t+\t2_167\t+\t*\n";  // bnd9

    // Additional events at the ends of chr1
    // del0 creates no edge
    // del4 creates no edge
    out << "S\tins0\tTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n";  // ins0
    out << "L\tins0\t+\t1_1\t+\t*\n";  // ins0
    out << "S\tins5\tTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n";  // ins5
    out << "L\t1_130\t+\tins5\t+\t*\n";  // ins5
    out << "L\t1_12\t+\t1_1\t+\t*\n";  // dup2
    out << "L\t1_130\t+\t1_130\t+\t*\n";  // dup3
    out << "L\t1_1\t-\t1_17\t+\t*\n";  // inv5
    out << "L\t1_129\t+\t1_130\t-\t*\n";  // inv6
    out << "S\trep2\tTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n";  // rep2
    out << "L\trep2\t+\t1_18\t+\t*\n";  // rep2
    out << "S\trep3\tTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n";  // rep3
    out << "L\t1_129\t+\trep3\t+\t*\n";  // rep3
    out << "L\t2_285\t+\t1_1\t+\t*\n";  // bnd21

    // Reference paths:
    //
    // chr1: 1_1+,1_2+,1_3+,1_4+,1_6+,1_7+,1_8+,1_9+,1_10+,1_11+,1_12+,1_13+,1_14+,1_17+,1_18+,1_19+,1_20+,1_21+,1_22+,1_23+,1_24+,1_25+,1_26+,1_28+,1_29+
    // chr1: 1_129+,1_130+
    //
    // chr2: 2_40+,2_50+,2_55+,2_58+,2_61+
    // chr2: 2_166+,2_167+,2_172+
    // chr2: 2_285+
    //
    // chr3: 3_60+,3_67+
    //
    // Paths covering each SV:
    if (print_paths_del) {
        out << "P\tdel1\t1_1+,1_2+,1_9+,1_10+,1_11+,1_12+,1_13+,1_14+,1_17+,1_18+,1_19+,1_20+,1_21+,1_22+,1_23+,1_24+,1_25+,1_26+,1_28+,1_29+\t*\n";
        out << "P\tdel2\t1_1+,1_2+,1_3+,1_4+,1_13+,1_14+,1_17+,1_18+,1_19+,1_20+,1_21+,1_22+,1_23+,1_24+,1_25+,1_26+,1_28+,1_29+\t*\n";
        //out << "P\tdel3\t1_1+,1_2+,1_3+,1_4+,1_6+,1_7+,1_17+,1_18+,1_19+,1_20+,1_21+,1_22+,1_23+,1_24+,1_25+,1_26+,1_28+,1_29+\t*\n";  // Commented out to test reverse orientation only
        out << "P\tdel3_rev\t1_29-,1_28-,1_26-,1_25-,1_24-,1_23-,1_22-,1_21-,1_20-,1_19-,1_18-,1_17-,1_7-,1_6-,1_4-,1_3-,1_2-,1_1-\t*\n";
        // No path covers del0
        // No path covers del4
    }
    if (print_paths_inv) {
        out << "P\tinv1\t1_1+,1_2+,1_3+,1_18-,1_17-,1_14-,1_13-,1_12-,1_11-,1_10-,1_9-,1_8-,1_7-,1_6-,1_4-,1_19+,1_20+,1_21+,1_22+,1_23+,1_24+,1_25+,1_26+,1_28+,1_29+\t*\n";
        out << "P\tinv2\t1_1+,1_2+,1_3+,1_4+,1_6+,1_21-,1_20-,1_19-,1_18-,1_17-,1_14-,1_13-,1_12-,1_11-,1_10-,1_9-,1_8-,1_7-,1_22+,1_23+,1_24+,1_25+,1_26+,1_28+,1_29+\t*\n";
        out << "P\tinv3\t1_1+,1_2+,1_3+,1_4+,1_6+,1_7+,1_8+,1_9+,1_10+,1_24-,1_23-,1_22-,1_21-,1_20-,1_19-,1_18-,1_17-,1_14-,1_13-,1_12-,1_11-,1_25+,1_26+,1_28+,1_29+\t*\n";
        //out << "P\tinv4\t1_1+,1_2+,1_3+,1_4+,1_6+,1_7+,1_8+,1_9+,1_10+,1_11+,1_12+,1_13+,1_14+,1_17+,1_18+,1_19+,1_20+,1_28-,1_26-,1_25-,1_24-,1_23-,1_22-,1_21-,1_29+\t*\n";  // Commented out to test reverse orientation only
        out << "P\tinv4_rev\t1_29-,1_21+,1_22+,1_23+,1_24+,1_25+,1_26+,1_28+,1_20-,1_19-,1_18-,1_17-,1_14-,1_13-,1_12-,1_11-,1_10-,1_9-,1_8-,1_7-,1_6-,1_4-,1_3-,1_2-,1_1-\t*\n";
        out << "P\tinv5\t1_14-,1_13-,1_12-,1_11-,1_10-,1_9-,1_8-,1_7-,1_6-,1_4-,1_3-,1_2-,1_1-,1_17+,1_18+,1_19+,1_20+,1_21+,1_22+,1_23+,1_24+,1_25+,1_26+,1_28+,1_29+\t*\n";
        out << "P\tinv6\t1_129+,1_130-\t*\n";
    }
    if (print_paths_ins) {
        out << "P\tins1\t1_1+,1_2+,1_3+,1_4+,ins1+,1_6+,1_7+,1_8+,1_9+,1_10+,1_11+,1_12+,1_13+,1_14+,1_17+,1_18+,1_19+,1_20+,1_21+,1_22+,1_23+,1_24+,1_25+,1_26+,1_28+,1_29+\t*\n";
        out << "P\tins2\t1_1+,1_2+,1_3+,1_4+,1_6+,1_7+,1_8+,1_9+,ins2+,1_10+,1_11+,1_12+,1_13+,1_14+,1_17+,1_18+,1_19+,1_20+,1_21+,1_22+,1_23+,1_24+,1_25+,1_26+,1_28+,1_29+\t*\n";
        out << "P\tins3\t1_1+,1_2+,1_3+,1_4+,1_6+,1_7+,1_8+,1_9+,1_10+,1_11+,ins3+,1_12+,1_13+,1_14+,1_17+,1_18+,1_19+,1_20+,1_21+,1_22+,1_23+,1_24+,1_25+,1_26+,1_28+,1_29+\t*\n";
        out << "P\tins3_rev\t1_29-,1_28-,1_26-,1_25-,1_24-,1_23-,1_22-,1_21-,1_20-,1_19-,1_18-,1_17-,1_14-,1_13-,1_12-,ins3-,1_11-,1_10-,1_9-,1_8-,1_7-,1_6-,1_4-,1_3-,1_2-,1_1-\t*\n";
        out << "P\tins0\tins0+,1_1+,1_2+,1_3+,1_4+,1_6+,1_7+,1_8+,1_9+,1_10+,1_11+,1_12+,1_13+,1_14+,1_17+,1_18+,1_19+,1_20+,1_21+,1_22+,1_23+,1_24+,1_25+,1_26+,1_28+,1_29+\t*\n";
        out << "P\tins5\t1_129+,1_130+,ins5+\t*\n";
    }
    if (print_paths_rep) {
        out << "P\trep1\t1_1+,1_2+,1_3+,1_4+,1_6+,1_7+,1_8+,rep1+,1_23+,1_24+,1_25+,1_26+,1_28+,1_29+\t*\n";
        out << "P\trep1_rev\t1_29-,1_28-,1_26-,1_25-,1_24-,1_23-,rep1-,1_8-,1_7-,1_6-,1_4-,1_3-,1_2-,1_1-\t*\n";
        out << "P\trep2\trep2+,1_18+,1_19+,1_20+,1_21+,1_22+,1_23+,1_24+,1_25+,1_26+,1_28+,1_29+\t*\n";
        out << "P\trep3\t1_129+,rep3+\t*\n";
    }
    if (print_paths_dup || print_paths_bnd) {
        out << "P\tdup1_bnd19\t1_1+,1_2+,1_3+,1_4+,1_6+,1_7+,1_8+,1_9+,1_10+,1_11+,1_12+,1_13+,1_14+,1_17+,1_18+,1_19+,1_20+,1_21+,1_22+,1_23+,1_24+,1_25+,1_26+,1_2+,1_3+,1_4+,1_6+,1_7+,1_8+,1_9+,1_10+,1_11+,1_12+,1_13+,1_14+,1_17+,1_18+,1_19+,1_20+,1_21+,1_22+,1_23+,1_24+,1_25+,1_26+,1_28+,1_29+\t*\n";
        out << "P\tdup1_bnd19_rev\t1_29-,1_28-,1_26-,1_25-,1_24-,1_23-,1_22-,1_21-,1_20-,1_19-,1_18-,1_17-,1_14-,1_13-,1_12-,1_11-,1_10-,1_9-,1_8-,1_7-,1_6-,1_4-,1_3-,1_2-,1_26-,1_25-,1_24-,1_23-,1_22-,1_21-,1_20-,1_19-,1_18-,1_17-,1_14-,1_13-,1_12-,1_11-,1_10-,1_9-,1_8-,1_7-,1_6-,1_4-,1_3-,1_2-,1_1-\t*\n";
    }
    if (print_paths_dup) {
        out << "P\tdup2\t1_1+,1_2+,1_3+,1_4+,1_6+,1_7+,1_8+,1_9+,1_10+,1_11+,1_12+,1_1+,1_2+,1_3+,1_4+,1_6+,1_7+,1_8+,1_9+,1_10+,1_11+,1_12+,1_13+,1_14+,1_17+,1_18+,1_19+,1_20+,1_21+,1_22+,1_23+,1_24+,1_25+,1_26+,1_28+,1_29+\t*\n";
        out << "P\tdup3\t1_129+,1_130+,1_130+\t*\n";
    }
    if (print_paths_bnd) {
        out << "P\tbnd1\t1_29-,1_28-,1_26-,1_25-,1_24-,1_23-,1_22-,1_21-,1_20-,1_19-,1_18-,1_17-,1_14-,1_13-,1_12-,1_11-,1_10-,1_9-,1_8-,1_7-,1_6-,1_4-,1_3-,2_40-\t*\n";
        out << "P\tbnd2\t1_1+,1_2+,1_3+,1_4+,1_6+,1_7+,2_50-,2_40-\t*\n";
        out << "P\tbnd3\t1_29-,1_28-,1_26-,1_25-,1_24-,1_23-,1_22-,1_21-,1_20-,1_19-,1_18-,1_17-,1_14-,1_13-,1_12-,1_11-,1_10-,1_9-,1_8-,2_58+,2_61+\t*\n";
        out << "P\tbnd4\t1_1+,1_2+,1_3+,1_4+,1_6+,1_7+,1_8+,1_9+,1_10+,1_11+,1_12+,1_13+,1_14+,1_17+,1_18+,1_19+,2_61+\t*\n";
        out << "P\tbnd5\t1_29-,1_28-,1_26-,1_25-,1_24-,1_23-,1_22-,1_21-,1_20-,bnd5-,2_172+ \t*\n";
        out << "P\tbnd6\t1_29-,1_28-,1_26-,1_25-,1_24-,bnd6-,3_67-,3_60-\t*\n";
        out << "P\tbnd7\t1_129+,bnd7+,3_67+\t*\n";
        out << "P\tbnd8\t1_129+,1_130+,bnd8+,2_285-\t*\n";
        out << "P\tbnd9\t1_1+,1_2+,1_3+,1_4+,1_6+,1_7+,1_8+,1_9+,1_10+,1_11+,1_12+,1_13+,2_167+,2_172+\t*\n";
        out << "P\tbnd10\tbnd10+,1_18+,1_19+,1_20+,1_21+,1_22+,1_23+,1_24+,1_25+,1_26+,1_28+,1_29+\t*\n";
        out << "P\tbnd11\t1_1+,1_2+,1_3+,1_4+,1_6+,1_7+,1_8+,1_9+,1_10+,1_11+,1_12+,1_13+,1_14+,1_17+,1_18+,1_19+,1_20+,1_21+,1_22+,1_23+,1_24+,1_25+,bnd11+\t*\n";
        out << "P\tbnd12\tbnd12+,1_26+,1_28+,1_29+\t*\n";
        out << "P\tbnd21\t1_29-,1_28-,1_26-,1_25-,1_24-,1_23-,1_22-,1_21-,1_20-,1_19-,1_18-,1_17-,1_14-,1_13-,1_12-,1_11-,1_10-,1_9-,1_8-,1_7-,1_6-,1_4-,1_3-,1_2-,1_1-,2_285-\t*\n";
    }
}


void print_truth_vcf_paths(ofstream& out) {
    const string QUAL = "60";  // Arbitrary
    const string FILTER = "PASS";
    const string FORMAT = "GT";
    const string SAMPLE = "0/1";  // Arbitrary
    const string SUFFIX_REPLACEMENT = QUAL+VcfReader::VCF_SEPARATOR+FILTER+VcfReader::VCF_SEPARATOR+VcfReader::VCF_MISSING_CHAR+VcfReader::VCF_SEPARATOR+FORMAT+VcfReader::VCF_SEPARATOR+SAMPLE+'\n';
    const string SUFFIX_BND = QUAL+VcfReader::VCF_SEPARATOR+FILTER+VcfReader::VCF_SEPARATOR+"SVTYPE=BND"+VcfReader::VCF_SEPARATOR+FORMAT+VcfReader::VCF_SEPARATOR+SAMPLE+'\n';

    out << "chr1\t1\tdel1\tAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tAAAAAAAAAAAAAAAAAAAAAA\t" << SUFFIX_REPLACEMENT;
    out << "chr1\t1\tdel2\tAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tAAAAAAAAAAAAAAAAAAAAA\t" << SUFFIX_REPLACEMENT;
    //out << "chr1\t1\tdel3\tAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tAAAAAAAAAAAAAAAAAAA\t" << SUFFIX_REPLACEMENT;  // Commented out to test reverse orientation only
    out << "chr1\t1\tdel3_rev\tAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tAAAAAAAAAAAAAAAAAAA\t" << SUFFIX_REPLACEMENT;
    out << "chr1\t1\tinv1\tAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tAAATTTTTTTTTTTTTTTAAAAAAAAAA\t" << SUFFIX_REPLACEMENT;
    out << "chr1\t1\tinv2\tAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tAAAAAATTTTTTTTTTTTTTTAAAAAAA\t" << SUFFIX_REPLACEMENT;
    out << "chr1\t1\tinv3\tAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tAAAAAAAAAATTTTTTTTTTTTTTAAAA\t" << SUFFIX_REPLACEMENT;
    //out << "chr1\t1\tinv4\tAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tAAAAAAAAAAAAAAAAAAAATTTTTTTT\t" << SUFFIX_REPLACEMENT;  // Commented out to test reverse orientation only
    out << "chr1\t1\tinv4_rev\tAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tAAAAAAAAAAAAAAAAAAAATTTTTTTT\t" << SUFFIX_REPLACEMENT;

    out << "chr1\t1\tins1\tAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tAAAAATTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAA\t" << SUFFIX_REPLACEMENT;
    out << "chr1\t1\tins2\tAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAA\t" << SUFFIX_REPLACEMENT;
    out << "chr1\t1\tins3\tAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAA\t" << SUFFIX_REPLACEMENT;
    out << "chr1\t1\tins3_rev\tAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAA\t" << SUFFIX_REPLACEMENT;
    out << "chr1\t1\trep1\tAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAA\t" << SUFFIX_REPLACEMENT;
    out << "chr1\t1\trep1_rev\tAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAA\t" << SUFFIX_REPLACEMENT;
    out << "chr1\t1\tdup1_bnd19\tAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t" << SUFFIX_REPLACEMENT;
    out << "chr1\t1\tdup1_bnd19_rev\tAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t" << SUFFIX_REPLACEMENT;

    out << "chr1\t29\tbnd1\tA\t]chr2:49]AAAAAAAAAAAAAAAAAAAAAAAAAAA\t" << SUFFIX_BND;
    out << "chr1\t1\tbnd2\tA\tAAAAAAAGGGGG]chr2:49]\t" << SUFFIX_BND;
    out << "chr1\t29\tbnd3\tA\t[chr2:61[GGGAAAAAAAAAAAAAAAAAAAAAA\t" << SUFFIX_BND;
    out << "chr1\t1\tbnd4\tA\tAAAAAAAAAAAAAAAAAAA[chr2:61[\t" << SUFFIX_BND;
    out << "chr1\t29\tbnd5\tA\t[chr2:172[TTTTTAAAAAAAAAA\t" << SUFFIX_BND;
    out << "chr1\t29\tbnd6\tA\t]chr3:59]GGGGGGGGTTTTTTTTTTAAAAAA\t" << SUFFIX_BND;
    out << "chr1\t129\tbnd7\tA\tATTTTTTTTTTTTTTT[chr3:60[\t" << SUFFIX_BND;
    out << "chr1\t129\tbnd8\tA\tAAATTTTTTTTTTTTTTTTTTTT]chr2:285]\t" << SUFFIX_BND;
    out << "chr1\t1\tbnd9\tA\tAAAAAAAAAAAAACCCCC[chr2:172[\t" << SUFFIX_BND;

    // Additional events at the ends of chr1
    // No path covers del0
    // No path covers del4
    out << "chr1\t29\tins0\tA\t.TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t" << SUFFIX_BND;
    out << "chr1\t129\tins5\tA\tAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT.\t" << SUFFIX_BND;
    out << "chr1\t29\tbnd10\tA\t.TTTTTAAAAAAAAAAAA\t" << SUFFIX_BND;
    out << "chr1\t29\tbnd12\tA\t.TTTTTTTTTTAAAA\t" << SUFFIX_BND;
    out << "chr1\t1\tbnd11\tA\tAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTT.\t" << SUFFIX_BND;
    out << "chr1\t1\tdup2\tAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t" << SUFFIX_REPLACEMENT;
    out << "chr1\t129\tdup3\tA\tAAA\t" << SUFFIX_REPLACEMENT;
    out << "chr1\t14\tinv5\tA\t[chr1:29[TTTTTTTTTTTTAAAAAAAAAAAAAA\t" << SUFFIX_BND;
    out << "chr1\t129\tinv6\tA\tA]chr1:131]\t" << SUFFIX_BND;
    out << "chr1\t29\trep2\tA\t.TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAA\t" << SUFFIX_BND;
    out << "chr1\t129\trep3\tA\tATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT.\t" << SUFFIX_BND;
    out << "chr1\t29\tbnd21\tA\t]chr2:285]AAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t" << SUFFIX_BND;
}


/**
 * Assumes that the symbolic INS is not in the list of VCF records.
 */
vector<pair<edge_t,uint32_t>> get_edge_record_map(const HashGraph& graph, const vector<string>& node_ids) {
    vector<pair<edge_t,uint32_t>> out;
    handle_t handle_from, handle_to;

    // dup1
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_26"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_2"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),0);  // dup1
    out.emplace_back(graph.edge_handle(handle_from,handle_to),1);  // bnd19
    out.emplace_back(graph.edge_handle(handle_from,handle_to),28);  // bnd20

    // del1
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_2"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_9"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),2);

    // inv1
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_18"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_3"))+1);
    handle_to=graph.flip(handle_to);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),3);
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_4"))+1);
    handle_from=graph.flip(handle_from);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_19"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),3);

    // bnd1
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"2_40"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_3"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),4);

    // del2
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_4"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_13"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),6);

    // ins1
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"ins1"))+1);
    handle_from=graph.flip(handle_from);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_4"))+1);
    handle_to=graph.flip(handle_to);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),7);
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"ins1"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_6"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),7);

    // inv2
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_6"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_21"))+1);
    handle_to=graph.flip(handle_to);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),9);
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_7"))+1);
    handle_from=graph.flip(handle_from);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_22"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),9);

    // del3
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_7"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_17"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),10);

    // bnd2
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_7"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"2_50"))+1);
    handle_to=graph.flip(handle_to);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),11);

    // rep1
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_8"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"rep1"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),12);
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"rep1"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_23"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),12);

    // bnd3
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_8"))+1);
    handle_from=graph.flip(handle_from);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"2_58"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),13);

    // ins2
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_9"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"ins2"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),14);
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"ins2"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_10"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),14);

    // inv3
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_10"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_24"))+1);
    handle_to=graph.flip(handle_to);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),15);
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_11"))+1);
    handle_from=graph.flip(handle_from);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_25"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),15);

    // ins3
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_11"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"ins3"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),16);
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"ins3"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_12"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),16);

    // bnd9
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_13"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"2_167"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),17);

    // bnd10
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"bnd10"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_18"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),20);

    // bnd4
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_19"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"2_61"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),21);

    // inv4
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_20"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_28"))+1);
    handle_to=graph.flip(handle_to);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),22);
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_21"))+1);
    handle_from=graph.flip(handle_from);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_29"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),22);

    // bnd5
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"bnd5"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_20"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),23);
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"2_172"))+1);
    handle_from=graph.flip(handle_from);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"bnd5"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),23);

    // bnd6
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"3_67"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"bnd6"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),25);
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"bnd6"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_24"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),25);

    // bnd11
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_25"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"bnd11"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),26);

    // bnd12
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"bnd12"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_26"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),27);

    // bnd7
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_129"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"bnd7"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),29);
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"bnd7"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"3_67"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),29);

    // bnd8
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_130"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"bnd8"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),30);  // Excluding symbolic INS
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"bnd8"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"2_285"))+1);
    handle_to=graph.flip(handle_to);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),30);  // Excluding symbolic INS

    // Additional events at the ends of chr1
    // No edge is associated with del0
    // No edge is associated with del4
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"ins0"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_1"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),34);
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_130"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"ins5"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),35);
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_12"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_1"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),36);
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_130"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_130"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),37);
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_1"))+1);
    handle_from=graph.flip(handle_from);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_17"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),38);
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_129"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_130"))+1);
    handle_to=graph.flip(handle_to);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),39);
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"rep2"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_18"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),40);
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_129"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"rep3"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),41);
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"2_285"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1_1"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),42);

    return out;
}


void node_to_chromosome_insert(VariantGraph& graph, const vector<string>& node_labels) {
    // chr1
    graph.node_to_chromosome_insert("1_1",node_labels,"chr1",0);
    graph.node_to_chromosome_insert("1_2",node_labels,"chr1",1);
    graph.node_to_chromosome_insert("1_3",node_labels,"chr1",2);
    graph.node_to_chromosome_insert("1_4",node_labels,"chr1",3);
    graph.node_to_chromosome_insert("1_6",node_labels,"chr1",5);
    graph.node_to_chromosome_insert("1_7",node_labels,"chr1",6);
    graph.node_to_chromosome_insert("1_8",node_labels,"chr1",7);
    graph.node_to_chromosome_insert("1_9",node_labels,"chr1",8);
    graph.node_to_chromosome_insert("1_10",node_labels,"chr1",9);
    graph.node_to_chromosome_insert("1_11",node_labels,"chr1",10);
    graph.node_to_chromosome_insert("1_12",node_labels,"chr1",11);
    graph.node_to_chromosome_insert("1_13",node_labels,"chr1",12);
    graph.node_to_chromosome_insert("1_14",node_labels,"chr1",13);
    graph.node_to_chromosome_insert("1_17",node_labels,"chr1",16);
    graph.node_to_chromosome_insert("1_18",node_labels,"chr1",17);
    graph.node_to_chromosome_insert("1_19",node_labels,"chr1",18);
    graph.node_to_chromosome_insert("1_20",node_labels,"chr1",19);
    graph.node_to_chromosome_insert("1_21",node_labels,"chr1",20);
    graph.node_to_chromosome_insert("1_22",node_labels,"chr1",21);
    graph.node_to_chromosome_insert("1_23",node_labels,"chr1",22);
    graph.node_to_chromosome_insert("1_24",node_labels,"chr1",23);
    graph.node_to_chromosome_insert("1_25",node_labels,"chr1",24);
    graph.node_to_chromosome_insert("1_26",node_labels,"chr1",25);
    graph.node_to_chromosome_insert("1_28",node_labels,"chr1",27);
    graph.node_to_chromosome_insert("1_29",node_labels,"chr1",28);

    graph.node_to_chromosome_insert("1_129",node_labels,"chr1",69);
    graph.node_to_chromosome_insert("1_130",node_labels,"chr1",129);

    // chr2
    graph.node_to_chromosome_insert("2_40",node_labels,"chr2",39);
    graph.node_to_chromosome_insert("2_50",node_labels,"chr2",49);
    graph.node_to_chromosome_insert("2_55",node_labels,"chr2",54);
    graph.node_to_chromosome_insert("2_58",node_labels,"chr2",57);
    graph.node_to_chromosome_insert("2_61",node_labels,"chr2",60);

    graph.node_to_chromosome_insert("2_166",node_labels,"chr2",131);
    graph.node_to_chromosome_insert("2_167",node_labels,"chr2",166);
    graph.node_to_chromosome_insert("2_172",node_labels,"chr2",171);

    graph.node_to_chromosome_insert("2_285",node_labels,"chr2",249);

    // chr3
    graph.node_to_chromosome_insert("3_60",node_labels,"chr3",39);
    graph.node_to_chromosome_insert("3_67",node_labels,"chr3",59);
}


unordered_map<string,string> get_chromosomes() {
    unordered_map<string,string> out;
    out.emplace("chr1","AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    out.emplace("chr2","CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
    out.emplace("chr3","GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
    return out;
}


/**
 * Zero-based, non-overlapping, [x..y).
 */
unordered_map<string,vector<interval_t>> get_tandem_track() {
    unordered_map<string,vector<interval_t>> out;

    vector<interval_t> vector1;
    vector1.emplace_back(25,29);
    vector1.emplace_back(30,35);
    vector1.emplace_back(37,41);
    vector1.emplace_back(79,97);
    vector1.emplace_back(99,130);
    out.insert({"chr1",vector1});

    vector<interval_t> vector2;
    vector2.emplace_back(58,62);
    vector2.emplace_back(69,80);
    vector2.emplace_back(141,150);
    vector2.emplace_back(155,165);
    vector2.emplace_back(171,176);
    vector2.emplace_back(179,192);
    vector2.emplace_back(259,272);
    vector2.emplace_back(279,285);
    out.insert({"chr2",vector2});

    vector<interval_t> vector3;
    vector3.emplace_back(49,67);
    out.insert({"chr3",vector3});

    return out;
}


int main(int argc, char* argv[]) {
    const path ROOT_DIR = path(argv[1]);

    const path INPUT_VCF = ROOT_DIR/"input.vcf";
    const path TRUTH_VCF = ROOT_DIR/"truth.vcf";
    const path TRUTH_GFA = ROOT_DIR/"truth.gfa";
    const path TEST_VCF = ROOT_DIR/"test.vcf";
    const path TEST_VCF_2 = ROOT_DIR/"test2.vcf";
    const path TEST_GFA = ROOT_DIR/"test.gfa";
    const uint16_t FLANK_LENGTH = 10;
    const uint16_t SIGNATURE_N_STEPS = 10;
    const uint16_t N_REUSE_TESTS = 10;

    ofstream input_vcf(INPUT_VCF.string());
    print_truth_vcf_header(input_vcf);
    print_truth_vcf(0,true/*Arbitrary*/,input_vcf);
    input_vcf.close();
    ofstream truth_gfa(TRUTH_GFA.string());
    print_truth_gfa(false,false,false,false,false,false,truth_gfa);
    truth_gfa.close();

    const unordered_map<string,string> chromosomes = get_chromosomes();
    const unordered_map<string,vector<interval_t>> tandem_track = get_tandem_track();
    string command;

    vector<VcfRecord> records;
    VcfReader reader(INPUT_VCF);
    reader.for_record_in_vcf([&](VcfRecord& record) {
        if ( (record.sv_type==VcfReader::TYPE_INSERTION && record.is_symbolic) ||
             ((record.sv_type==VcfReader::TYPE_DELETION || record.sv_type==VcfReader::TYPE_INVERSION || record.sv_type==VcfReader::TYPE_DUPLICATION || record.sv_type==VcfReader::TYPE_REPLACEMENT) && record.sv_length==UINT32_MAX)
           ) return;
        records.push_back(record);
    });
    const vector<string> caller_ids = get_caller_ids();
    auto rd = std::random_device {};
    VariantGraph graph(chromosomes,tandem_track);
    for (uint16_t i=0; i<N_REUSE_TESTS; i++) {  // The same VariantGraph class can be reused multiple times
        cerr << "----- REUSE TEST " << std::to_string((i+1)) << " ------\n";
        auto rng = std::default_random_engine { rd() };
        vector<VcfRecord> records_prime = records;
        std::shuffle(std::begin(records_prime),std::end(records_prime),rng);
        graph.build(records_prime,FLANK_LENGTH,false,caller_ids);

        cerr << "Testing print_supported_vcf_records() (all records)...\n";
        ofstream truth_vcf(TRUTH_VCF.string());
        print_truth_vcf_header(truth_vcf);
        print_truth_vcf(1,true/*Arbitrary*/,truth_vcf);
        truth_vcf.close();
        ofstream test_vcf(TEST_VCF.string());
        print_truth_vcf_header(test_vcf);
        graph.print_supported_vcf_records(test_vcf,true,caller_ids);
        test_vcf.close();
        command.clear(); command.append("bcftools view --no-header "+TRUTH_VCF.string()+" | sort > tmp1.txt"); run_command(command);
        command.clear(); command.append("bcftools view --no-header "+TEST_VCF.string()+" | sort > tmp2.txt"); run_command(command);
        command.clear(); command.append("diff --brief tmp1.txt tmp2.txt"); run_command(command);

        cerr << "Testing to_gfa(): node sequences...\n";
        graph.to_gfa(TEST_GFA);
        command.clear(); command.append("grep ^S "+TRUTH_GFA.string()+" | cut -f 1,3 | sort > tmp1.txt"); run_command(command);
        command.clear(); command.append("grep ^S "+TEST_GFA.string()+" | cut -f 1,3 | sort > tmp2.txt"); run_command(command);
        command.clear(); command.append("diff --brief tmp1.txt tmp2.txt"); run_command(command);

        cerr << "Testing to_gfa(): n. edges...\n";
        command.clear(); command.append("grep ^L "+TRUTH_GFA.string()+" | wc -l > tmp1.txt"); run_command(command);
        command.clear(); command.append("grep ^L "+TEST_GFA.string()+" | wc -l > tmp2.txt"); run_command(command);
        command.clear(); command.append("diff --brief tmp1.txt tmp2.txt"); run_command(command);

        cerr << "Testing to_gfa(): local topology (" << SIGNATURE_N_STEPS << " steps)...\n";
        graph.print_graph_signature(SIGNATURE_N_STEPS,"tmp1.txt");
        graph.load_gfa(TRUTH_GFA);
        graph.print_graph_signature(SIGNATURE_N_STEPS,"tmp2.txt");
        command.clear(); command.append("diff --brief tmp1.txt tmp2.txt"); run_command(command);
    }

    cerr << "Testing print_supported_vcf_records() (paths)...\n";
    ofstream truth_vcf(TRUTH_VCF.string());
    print_truth_vcf_header(truth_vcf);
    print_truth_vcf(2,false/*Arbitrary*/,truth_vcf);
    truth_vcf.close();
    graph.build(records,FLANK_LENGTH,false,caller_ids);

    cerr << "DEL...\n";
    truth_gfa.clear(); truth_gfa.open(TRUTH_GFA.string());
    print_truth_gfa(true,false,false,false,false,false,truth_gfa);
    truth_gfa.close();
    vector<string> node_labels = graph.load_gfa(TRUTH_GFA.string());
    vector<pair<edge_t,uint32_t>> edge_record_map = get_edge_record_map(graph.graph,node_labels);
    graph.load_edge_record_map(edge_record_map,get_n_vcf_records());
    ofstream test_vcf(TEST_VCF.string());
    print_truth_vcf_header(test_vcf);
    graph.print_supported_vcf_records(test_vcf,false,caller_ids);
    test_vcf.close();
    command.clear(); command.append("bcftools view --no-header "+TRUTH_VCF.string()+" | grep 'SVTYPE=DEL' | sort > tmp1.txt"); run_command(command);
    command.clear(); command.append("bcftools view --no-header "+TEST_VCF.string()+" | sort > tmp2.txt"); run_command(command);
    command.clear(); command.append("diff --brief tmp1.txt tmp2.txt"); run_command(command);

    cerr << "INV...\n";
    truth_gfa.clear(); truth_gfa.open(TRUTH_GFA.string());
    print_truth_gfa(false,true,false,false,false,false,truth_gfa);
    truth_gfa.close();
    node_labels=graph.load_gfa(TRUTH_GFA.string());
    test_vcf.clear(); test_vcf.open(TEST_VCF.string());
    print_truth_vcf_header(test_vcf);
    graph.print_supported_vcf_records(test_vcf,false,caller_ids);
    test_vcf.close();
    command.clear(); command.append("bcftools view --no-header "+TRUTH_VCF.string()+" | grep 'SVTYPE=INV' | sort > tmp1.txt"); run_command(command);
    command.clear(); command.append("bcftools view --no-header "+TEST_VCF.string()+" | sort > tmp2.txt"); run_command(command);
    command.clear(); command.append("diff --brief tmp1.txt tmp2.txt"); run_command(command);

    cerr << "INS...\n";
    truth_gfa.clear(); truth_gfa.open(TRUTH_GFA.string());
    print_truth_gfa(false,false,true,false,false,false,truth_gfa);
    truth_gfa.close();
    node_labels=graph.load_gfa(TRUTH_GFA.string());
    test_vcf.clear(); test_vcf.open(TEST_VCF.string());
    print_truth_vcf_header(test_vcf);
    graph.print_supported_vcf_records(test_vcf,false,caller_ids);
    test_vcf.close();
    command.clear(); command.append("bcftools view --no-header "+TRUTH_VCF.string()+" | grep 'SVTYPE=INS' | sort > tmp1.txt"); run_command(command);
    command.clear(); command.append("bcftools view --no-header "+TEST_VCF.string()+" | sort > tmp2.txt"); run_command(command);
    command.clear(); command.append("diff --brief tmp1.txt tmp2.txt"); run_command(command);

    cerr << "REP...\n";
    truth_gfa.clear(); truth_gfa.open(TRUTH_GFA.string());
    print_truth_gfa(false,false,false,true,false,false,truth_gfa);
    truth_gfa.close();
    node_labels=graph.load_gfa(TRUTH_GFA.string());
    test_vcf.clear(); test_vcf.open(TEST_VCF.string());
    print_truth_vcf_header(test_vcf);
    graph.print_supported_vcf_records(test_vcf,false,caller_ids);
    test_vcf.close();
    command.clear(); command.append("bcftools view --no-header "+TRUTH_VCF.string()+" | grep rep | sort > tmp1.txt"); run_command(command);
    command.clear(); command.append("bcftools view --no-header "+TEST_VCF.string()+" | sort > tmp2.txt"); run_command(command);
    command.clear(); command.append("diff --brief tmp1.txt tmp2.txt"); run_command(command);

    cerr << "DUP...\n";
    truth_gfa.clear(); truth_gfa.open(TRUTH_GFA.string());
    print_truth_gfa(false,false,false,false,true,false,truth_gfa);
    truth_gfa.close();
    node_labels=graph.load_gfa(TRUTH_GFA.string());
    test_vcf.clear(); test_vcf.open(TEST_VCF.string());
    print_truth_vcf_header(test_vcf);
    graph.print_supported_vcf_records(test_vcf,false,caller_ids);
    test_vcf.close();
    command.clear(); command.append("bcftools view --no-header "+TRUTH_VCF.string()+" | grep -E 'SVTYPE=DUP|bnd19|bnd20' | sort | cut -f 1,2,3,4,5,6,7,9,10 > tmp1.txt"); run_command(command);
    command.clear(); command.append("bcftools view --no-header "+TEST_VCF.string()+" | sort | cut -f 1,2,3,4,5,6,7,9,10 > tmp2.txt"); run_command(command);
    command.clear(); command.append("diff --brief tmp1.txt tmp2.txt"); run_command(command);

    cerr << "BND...\n";
    truth_gfa.clear(); truth_gfa.open(TRUTH_GFA.string());
    print_truth_gfa(false,false,false,false,false,true,truth_gfa);
    truth_gfa.close();
    node_labels=graph.load_gfa(TRUTH_GFA.string());
    test_vcf.clear(); test_vcf.open(TEST_VCF.string());
    print_truth_vcf_header(test_vcf);
    graph.print_supported_vcf_records(test_vcf,false,caller_ids);
    test_vcf.close();
    command.clear(); command.append("bcftools view --no-header "+TRUTH_VCF.string()+" | grep -E 'SVTYPE=BND|dup1' | sort | cut -f 1,2,3,4,5,6,7,9,10 > tmp1.txt"); run_command(command);
    command.clear(); command.append("bcftools view --no-header "+TEST_VCF.string()+" | sort | cut -f 1,2,3,4,5,6,7,9,10 > tmp2.txt"); run_command(command);
    command.clear(); command.append("diff --brief tmp1.txt tmp2.txt"); run_command(command);

    cerr << "All types...\n";
    truth_gfa.clear(); truth_gfa.open(TRUTH_GFA.string());
    print_truth_gfa(true,true,true,true,true,true,truth_gfa);
    truth_gfa.close();
    node_labels=graph.load_gfa(TRUTH_GFA.string());
    test_vcf.clear(); test_vcf.open(TEST_VCF.string());
    print_truth_vcf_header(test_vcf);
    graph.print_supported_vcf_records(test_vcf,false,caller_ids);
    test_vcf.close();
    command.clear(); command.append("bcftools view --no-header "+TRUTH_VCF.string()+" | sort | cut -f 1,2,3,4,5,6,7,9,10 > tmp1.txt"); run_command(command);
    command.clear(); command.append("bcftools view --no-header "+TEST_VCF.string()+" | sort | cut -f 1,2,3,4,5,6,7,9,10 > tmp2.txt"); run_command(command);
    command.clear(); command.append("diff --brief tmp1.txt tmp2.txt"); run_command(command);

    cerr << "Testing paths_to_vcf_records()...\n";
    test_vcf.clear(); test_vcf.open(TEST_VCF.string());
    print_truth_vcf_header(test_vcf);
    graph.node_to_chromosome_clear(); node_to_chromosome_insert(graph,node_labels);
    graph.paths_to_vcf_records(test_vcf);
    test_vcf.close();
    truth_vcf.clear(); truth_vcf.open(TRUTH_VCF.string());
    print_truth_vcf_header(truth_vcf);
    print_truth_vcf_paths(truth_vcf);
    truth_vcf.close();
    command.clear(); command.append("bcftools view --no-header "+TRUTH_VCF.string()+" | cut -f 1,2,4,5 | sort > tmp1.txt"); run_command(command);
    command.clear(); command.append("bcftools view --no-header "+TEST_VCF.string()+" | cut -f 1,2,4,5 | sort > tmp2.txt"); run_command(command);
    command.clear(); command.append("diff --brief tmp1.txt tmp2.txt"); run_command(command);

    cerr << "Loading back the output of paths_to_vcf_records()...\n";
    records.clear();
    VcfReader reader2(TEST_VCF);
    reader2.for_record_in_vcf([&](VcfRecord& record) {
        if ( (record.sv_type==VcfReader::TYPE_INSERTION && record.is_symbolic) ||
             ((record.sv_type==VcfReader::TYPE_DELETION || record.sv_type==VcfReader::TYPE_INVERSION || record.sv_type==VcfReader::TYPE_DUPLICATION || record.sv_type==VcfReader::TYPE_REPLACEMENT) && record.sv_length==UINT32_MAX)
           ) return;
        records.push_back(record);
    });
    graph.build(records,FLANK_LENGTH,false,caller_ids);
    ofstream test_vcf2(TEST_VCF_2.string());
    print_truth_vcf_header(test_vcf2);
    graph.print_supported_vcf_records(test_vcf2,true,caller_ids);
    test_vcf2.close();
    command.clear(); command.append("bcftools view --no-header "+TEST_VCF.string()+" | sort > tmp1.txt"); run_command(command);
    command.clear(); command.append("bcftools view --no-header "+TEST_VCF_2.string()+" | sort > tmp2.txt"); run_command(command);
    command.clear(); command.append("diff --brief tmp1.txt tmp2.txt"); run_command(command);

    cerr << "Removing temporary files...\n";
    command.clear(); command.append("rm -f truth*.vcf test*.vcf truth*.gfa test*.gfa tmp*.txt"); run_command(command);
}