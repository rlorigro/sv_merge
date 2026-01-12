#include "VcfReader.hpp"
#include "VariantGraph.hpp"
#include "misc.hpp"

using std::stoi;
using sv_merge::interval_t;
using sv_merge::run_command;
using sv_merge::VcfReader;
using sv_merge::VariantGraph;
using bdsg::step_handle_t;
using bdsg::HashGraph;

#include <iostream>
#include <algorithm>
#include <random>


using std::ofstream;


void print_truth_vcf_header(ofstream& out) {
    out << "##fileformat=VCFv4.2\n";
    out << "##contig=<ID=chr1,length=131>\n";
    out << R"(##ALT=<ID=DUP,Description="Duplication">)" << '\n';
    out << R"(##ALT=<ID=CNV,Description="Copy-number variant">)" << '\n';
    out << R"(##ALT=<ID=INV,Description="Inversion">)" << '\n';
    out << R"(##FILTER=<ID=PASS,Description="All filters passed">)" << '\n';
    out << R"(##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">)" << '\n';
    out << R"(##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variation">)" << '\n';
    out << R"(##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variation">)" << '\n';
    out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
}


/**
 * Remark: BNDs are not tested; we leave this to the future, since BND handling should be tested more thoroughly in the
 * entire codebase.
 */
void print_truth_vcf(size_t dataset, ofstream& out) {
    const string QUAL = "60";  // Arbitrary
    const string FILTER = "PASS";
    const string FORMAT = "GT";
    const string GT = "0/1";  // Arbitrary
    const string INFIX = QUAL+"\t"+FILTER;
    const string SUFFIX = FORMAT+"\t"+GT+"\n";

    if (dataset==0) {
        // INS section
        out << "chr1\t10\tdel1\tAAAAAA\tA\t" << INFIX << "\tSVTYPE=DEL;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t15\tdel2\tAAAAAA\tA\t" << INFIX << "\tSVTYPE=DEL;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t15\tins1\tA\tACCCCC\t" << INFIX << "\tSVTYPE=INS;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t15\tins2\tA\tAGGGGG\t" << INFIX << "\tSVTYPE=INS;SVLEN=5;\t" << SUFFIX;

        out << "chr1\t31\trep1\tCCCCC\tAAAAA\t" << INFIX << "\tSVLEN=5;\t" << SUFFIX;
        out << "chr1\t31\trep2\tCCCCC\tGGGGG\t" << INFIX << "\tSVLEN=5;\t" << SUFFIX;
        out << "chr1\t35\tins3\tC\tCAAAAA\t" << INFIX << "\tSVTYPE=INS;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t35\tins4\tC\tCGGGGG\t" << INFIX << "\tSVTYPE=INS;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t36\trep3\tCCCCC\tAAAAA\t" << INFIX << "\tSVLEN=5;\t" << SUFFIX;
        out << "chr1\t36\trep4\tCCCCC\tGGGGG\t" << INFIX << "\tSVLEN=5;\t" << SUFFIX;

        out << "chr1\t50\tdup1\tG\t<DUP>\t" << INFIX << "\tSVTYPE=DUP;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t55\tins5\tG\tGAAAAA\t" << INFIX << "\tSVTYPE=INS;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t55\tins6\tG\tGCCCCC\t" << INFIX << "\tSVTYPE=INS;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t55\tdup2\tG\t<DUP>\t" << INFIX << "\tSVTYPE=DUP;SVLEN=5;\t" << SUFFIX;

        out << "chr1\t70\tinv1\tT\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t75\tins7\tT\tTAAAAA\t" << INFIX << "\tSVTYPE=INS;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t75\tins8\tT\tTGGGGG\t" << INFIX << "\tSVTYPE=INS;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t75\tinv2\tT\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=5;\t" << SUFFIX;

        // DEL section
        out << "chr1\t90\tdel3\tAAAAAA\tA\t" << INFIX << "\tSVTYPE=DEL;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t95\tdel4\tAAAAAA\tA\t" << INFIX << "\tSVTYPE=DEL;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t100\tdel5\tAAAAAA\tA\t" << INFIX << "\tSVTYPE=DEL;SVLEN=5;\t" << SUFFIX;

        out << "chr1\t116\trep5\tCCCCC\tAAAAA\t" << INFIX << "\tSVLEN=5;\t" << SUFFIX;
        out << "chr1\t116\trep6\tCCCCC\tGGGGG\t" << INFIX << "\tSVLEN=5;\t" << SUFFIX;
        out << "chr1\t120\tdel6\tCCCCCC\tC\t" << INFIX << "\tSVTYPE=DEL;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t126\trep7\tCCCCC\tAAAAA\t" << INFIX << "\tSVLEN=5;\t" << SUFFIX;
        out << "chr1\t126\trep8\tCCCCC\tGGGGG\t" << INFIX << "\tSVLEN=5;\t" << SUFFIX;

        out << "chr1\t140\tdup3\tG\t<DUP>\t" << INFIX << "\tSVTYPE=DUP;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t145\tdel7\tGGGGGG\tG\t" << INFIX << "\tSVTYPE=DEL;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t150\tdup4\tG\t<DUP>\t" << INFIX << "\tSVTYPE=DUP;SVLEN=5;\t" << SUFFIX;

        out << "chr1\t165\tinv3\tT\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t170\tdel8\tTTTTTT\tT\t" << INFIX << "\tSVTYPE=DEL;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t175\tinv4\tT\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=5;\t" << SUFFIX;

        // REP section
        out << "chr1\t190\tinv5\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t196\trep9\tAAAAA\tCCCCC\t" << INFIX << "\tSVLEN=5;\t" << SUFFIX;
        out << "chr1\t196\trep10\tAAAAA\tGGGGG\t" << INFIX << "\tSVLEN=5;\t" << SUFFIX;
        out << "chr1\t200\tinv6\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=5;\t" << SUFFIX;

        out << "chr1\t215\tdup5\tC\t<DUP>\t" << INFIX << "\tSVTYPE=DUP;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t221\trep11\tCCCCC\tAAAAA\t" << INFIX << "\tSVLEN=5;\t" << SUFFIX;
        out << "chr1\t221\trep12\tCCCCC\tGGGGG\t" << INFIX << "\tSVLEN=5;\t" << SUFFIX;
        out << "chr1\t225\tdup6\tC\t<DUP>\t" << INFIX << "\tSVTYPE=DUP;SVLEN=5;\t" << SUFFIX;

        // DUP section
        out << "chr1\t240\tinv7\tG\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t245\tdup7\tG\t<DUP>\t" << INFIX << "\tSVTYPE=DUP;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t250\tinv8\tG\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=5;\t" << SUFFIX;

        out << "chr1\t265\tdup8\tT\t<DUP>\t" << INFIX << "\tSVTYPE=DUP;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t270\tdup9\tT\t<DUP>\t" << INFIX << "\tSVTYPE=DUP;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t275\tdup9\tT\t<DUP>\t" << INFIX << "\tSVTYPE=DUP;SVLEN=5;\t" << SUFFIX;

        // INV section
        out << "chr1\t290\tinv9\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t295\tinv10\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t300\tinv11\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=5;\t" << SUFFIX;
    }
    else if (dataset==1) {
        // REP section
        out << "chr1\t16\trep1\tCCCCC\tAAAAA\t" << INFIX << "\tSVLEN=5;\t" << SUFFIX;
        out << "chr1\t16\trep2\tCCCCC\tGGGGG\t" << INFIX << "\tSVLEN=5;\t" << SUFFIX;
        out << "chr1\t21\trep3\tCCCCC\tAAAAA\t" << INFIX << "\tSVLEN=5;\t" << SUFFIX;
        out << "chr1\t21\trep4\tCCCCC\tGGGGG\t" << INFIX << "\tSVLEN=5;\t" << SUFFIX;
        out << "chr1\t26\trep5\tCCCCC\tAAAAA\t" << INFIX << "\tSVLEN=5;\t" << SUFFIX;
        out << "chr1\t26\trep6\tCCCCC\tGGGGG\t" << INFIX << "\tSVLEN=5;\t" << SUFFIX;

        // DEL section
        out << "chr1\t50\tdel1\tGGGGGG\tG\t" << INFIX << "\tSVTYPE=DEL;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t50\tins1\tG\tGCCCCC\t" << INFIX << "\tSVTYPE=INS;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t50\tins2\tG\tGTTTTT\t" << INFIX << "\tSVTYPE=INS;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t55\tins3\tG\tGCCCCC\t" << INFIX << "\tSVTYPE=INS;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t55\tins4\tG\tGTTTTT\t" << INFIX << "\tSVTYPE=INS;SVLEN=5;\t" << SUFFIX;

        // REP section
        out << "chr1\t65\tins5\tT\tTAAAAA\t" << INFIX << "\tSVTYPE=INS;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t65\tins6\tT\tTCCCCC\t" << INFIX << "\tSVTYPE=INS;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t66\trep7\tTTTTT\tAAAAA\t" << INFIX << "\tSVLEN=5;\t" << SUFFIX;
        out << "chr1\t66\trep8\tTTTTT\tCCCCC\t" << INFIX << "\tSVLEN=5;\t" << SUFFIX;
        out << "chr1\t70\tins7\tT\tTAAAAA\t" << INFIX << "\tSVTYPE=INS;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t70\tins8\tT\tTCCCCC\t" << INFIX << "\tSVTYPE=INS;SVLEN=5;\t" << SUFFIX;

        out << "chr1\t80\tdel2\tAAAAAA\tA\t" << INFIX << "\tSVTYPE=DEL;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t86\trep9\tAAAAA\tCCCCC\t" << INFIX << "\tSVLEN=5;\t" << SUFFIX;
        out << "chr1\t86\trep10\tAAAAA\tGGGGG\t" << INFIX << "\tSVLEN=5;\t" << SUFFIX;
        out << "chr1\t90\tdel3\tAAAAAA\tA\t" << INFIX << "\tSVTYPE=DEL;SVLEN=5;\t" << SUFFIX;

        // INV section
        out << "chr1\t105\tins9\tC\tCGGGGG\t" << INFIX << "\tSVTYPE=INS;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t105\tins10\tC\tCTTTTT\t" << INFIX << "\tSVTYPE=INS;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t105\tinv12\tC\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t110\tins11\tC\tCGGGGG\t" << INFIX << "\tSVTYPE=INS;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t110\tins12\tC\tCTTTTT\t" << INFIX << "\tSVTYPE=INS;SVLEN=5;\t" << SUFFIX;

        out << "chr1\t120\tdel4\tGGGGGG\tG\t" << INFIX << "\tSVTYPE=DEL;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t125\tinv13\tG\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t130\tdel5\tGGGGGG\tG\t" << INFIX << "\tSVTYPE=DEL;SVLEN=5;\t" << SUFFIX;

        out << "chr1\t146\trep11\tTTTTT\tAAAAA\t" << INFIX << "\tSVLEN=5;\t" << SUFFIX;
        out << "chr1\t146\trep12\tTTTTT\tCCCCC\t" << INFIX << "\tSVLEN=5;\t" << SUFFIX;
        out << "chr1\t150\tinv14\tT\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t156\trep13\tTTTTT\tAAAAA\t" << INFIX << "\tSVLEN=5;\t" << SUFFIX;
        out << "chr1\t156\trep14\tTTTTT\tCCCCC\t" << INFIX << "\tSVLEN=5;\t" << SUFFIX;

        out << "chr1\t170\tdup1\tA\t<DUP>\t" << INFIX << "\tSVTYPE=DUP;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t175\tinv15\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=5;\t" << SUFFIX;
        out << "chr1\t180\tdup2\tA\t<DUP>\t" << INFIX << "\tSVTYPE=DUP;SVLEN=5;\t" << SUFFIX;
    }
}


void print_truth_gfa(size_t dataset, ofstream& out, bool closure1, bool closure2) {
    if (dataset==0) {
        // INS section
        // Block 1
        out << "S\t1\tAAAAAAAAAA\n";
        out << "S\t2\tAAAAA\n";
        out << "S\tins1\tCCCCC\n";
        out << "S\tins2\tGGGGG\n";
        out << "S\t3\tAAAAA\n";
        out << "S\t4\tAAAAACCCCC\n";
        // Before closure
        out << "L\t1\t+\t2\t+\t*\n";
        out << "L\t2\t+\t3\t+\t*\n";
        out << "L\t3\t+\t4\t+\t*\n";
        out << "L\t1\t+\t3\t+\t*\n";
        out << "L\t2\t+\t4\t+\t*\n";
        out << "L\t2\t+\tins1\t+\t*\n";
        out << "L\t2\t+\tins2\t+\t*\n";
        out << "L\tins1\t+\t3\t+\t*\n";
        out << "L\tins2\t+\t3\t+\t*\n";
        if (closure1) {
            out << "L\t1\t+\t4\t+\t*\n";
            out << "L\t1\t+\tins1\t+\t*\n";
            out << "L\t1\t+\tins2\t+\t*\n";
            out << "L\tins1\t+\t4\t+\t*\n";
            out << "L\tins2\t+\t4\t+\t*\n";
            if (closure2) {
                // NOP
            }
        }

        // Block 2
        out << "S\t5\tCCCCC\n";
        out << "S\t6\tCCCCC\n";
        out << "S\t7\tCCCCCGGGGG\n";
        out << "S\trep1\tAAAAA\n";
        out << "S\trep2\tGGGGG\n";
        out << "S\tins3\tAAAAA\n";
        out << "S\tins4\tGGGGG\n";
        out << "S\trep3\tAAAAA\n";
        out << "S\trep4\tGGGGG\n";
        // Before closure
        out << "L\t4\t+\t5\t+\t*\n";
        out << "L\t5\t+\t6\t+\t*\n";
        out << "L\t6\t+\t7\t+\t*\n";
        out << "L\t4\t+\trep1\t+\t*\n";
        out << "L\t4\t+\trep2\t+\t*\n";
        out << "L\trep1\t+\t6\t+\t*\n";
        out << "L\trep2\t+\t6\t+\t*\n";
        out << "L\t5\t+\tins3\t+\t*\n";
        out << "L\t5\t+\tins4\t+\t*\n";
        out << "L\tins3\t+\t6\t+\t*\n";
        out << "L\tins4\t+\t6\t+\t*\n";
        out << "L\t5\t+\trep3\t+\t*\n";
        out << "L\t5\t+\trep4\t+\t*\n";
        out << "L\trep3\t+\t7\t+\t*\n";
        out << "L\trep4\t+\t7\t+\t*\n";
        if (closure1) {
            out << "L\trep1\t+\tins3\t+\t*\n";
            out << "L\trep1\t+\tins4\t+\t*\n";
            out << "L\trep1\t+\trep3\t+\t*\n";
            out << "L\trep1\t+\trep4\t+\t*\n";
            out << "L\trep2\t+\tins3\t+\t*\n";
            out << "L\trep2\t+\tins4\t+\t*\n";
            out << "L\trep2\t+\trep3\t+\t*\n";
            out << "L\trep2\t+\trep4\t+\t*\n";
            out << "L\tins3\t+\trep3\t+\t*\n";
            out << "L\tins3\t+\trep4\t+\t*\n";
            out << "L\tins4\t+\trep3\t+\t*\n";
            out << "L\tins4\t+\trep4\t+\t*\n";
            if (closure2) {
                // NOP
            }
        }

        // Block 3
        out << "S\t8\tGGGGG\n";
        out << "S\t9\tGGGGG\n";
        out << "S\t10\tGGGGGTTTTT\n";
        out << "S\tins5\tAAAAA\n";
        out << "S\tins6\tCCCCC\n";
        // Before closure
        out << "L\t7\t+\t8\t+\t*\n";
        out << "L\t8\t+\t9\t+\t*\n";
        out << "L\t9\t+\t10\t+\t*\n";
        out << "L\t8\t+\t8\t+\t*\n";
        out << "L\t8\t+\tins5\t+\t*\n";
        out << "L\t8\t+\tins6\t+\t*\n";
        out << "L\tins5\t+\t9\t+\t*\n";
        out << "L\tins6\t+\t9\t+\t*\n";
        out << "L\t9\t+\t9\t+\t*\n";
        if (closure1) {
            // NOP
        }

        // Block 4
        out << "S\t11\tTTTTT\n";
        out << "S\t12\tTTTTT\n";
        out << "S\t13\tTTTTTAAAAA\n";
        out << "S\tins7\tAAAAA\n";
        out << "S\tins8\tGGGGG\n";
        // Before closure
        out << "L\t10\t+\t11\t+\t*\n";
        out << "L\t11\t+\t12\t+\t*\n";
        out << "L\t12\t+\t13\t+\t*\n";
        out << "L\t10\t+\t11\t-\t*\n";
        out << "L\t11\t-\t12\t+\t*\n";
        out << "L\t11\t+\tins7\t+\t*\n";
        out << "L\t11\t+\tins8\t+\t*\n";
        out << "L\tins7\t+\t12\t+\t*\n";
        out << "L\tins8\t+\t12\t+\t*\n";
        out << "L\t11\t+\t12\t-\t*\n";
        out << "L\t12\t-\t13\t+\t*\n";
        if (closure1) {
            out << "L\t10\t+\tins7\t-\t*\n";
            out << "L\t10\t+\tins8\t-\t*\n";
            out << "L\t10\t+\t13\t+\t*\n";
            out << "L\t11\t-\tins7\t+\t*\n";
            out << "L\t11\t-\tins8\t+\t*\n";
            out << "L\t11\t-\t12\t-\t*\n";
            out << "L\tins7\t-\t13\t+\t*\n";
            out << "L\tins8\t-\t13\t+\t*\n";
            out << "L\tins7\t+\t12\t-\t*\n";
            out << "L\tins8\t+\t12\t-\t*\n";
            if (closure2) {
                // NOP
            }
        }

        // DEL section
        // Block 1
        out << "S\t14\tAAAAA\n";
        out << "S\t15\tAAAAA\n";
        out << "S\t16\tAAAAA\n";
        out << "S\t17\tAAAAACCCCC\n";
        // Before closure
        out << "L\t13\t+\t14\t+\t*\n";
        out << "L\t14\t+\t15\t+\t*\n";
        out << "L\t15\t+\t16\t+\t*\n";
        out << "L\t16\t+\t17\t+\t*\n";
        out << "L\t13\t+\t15\t+\t*\n";
        out << "L\t14\t+\t16\t+\t*\n";
        out << "L\t15\t+\t17\t+\t*\n";
        if (closure1) {
            out << "L\t13\t+\t16\t+\t*\n";
            out << "L\t14\t+\t17\t+\t*\n";
            if (closure2) {
                out << "L\t13\t+\t17\t+\t*\n";
            }
        }

        // Block 2
        out << "S\t18\tCCCCC\n";
        out << "S\trep5\tAAAAA\n";
        out << "S\trep6\tGGGGG\n";
        out << "S\t19\tCCCCC\n";
        out << "S\t20\tCCCCC\n";
        out << "S\trep7\tAAAAA\n";
        out << "S\trep8\tGGGGG\n";
        out << "S\t21\tCCCCCGGGGG\n";
        // Before closure
        out << "L\t17\t+\t18\t+\t*\n";
        out << "L\t18\t+\t19\t+\t*\n";
        out << "L\t19\t+\t20\t+\t*\n";
        out << "L\t20\t+\t21\t+\t*\n";
        out << "L\t17\t+\trep5\t+\t*\n";
        out << "L\t17\t+\trep6\t+\t*\n";
        out << "L\trep5\t+\t19\t+\t*\n";
        out << "L\trep6\t+\t19\t+\t*\n";
        out << "L\t18\t+\t20\t+\t*\n";
        out << "L\t19\t+\trep7\t+\t*\n";
        out << "L\t19\t+\trep8\t+\t*\n";
        out << "L\trep7\t+\t21\t+\t*\n";
        out << "L\trep8\t+\t21\t+\t*\n";
        if (closure1) {
            out << "L\trep5\t+\t20\t+\t*\n";
            out << "L\trep6\t+\t20\t+\t*\n";
            out << "L\t18\t+\trep7\t+\t*\n";
            out << "L\t18\t+\trep8\t+\t*\n";
            if (closure2) {
                out << "L\trep5\t+\trep7\t+\t*\n";
                out << "L\trep5\t+\trep8\t+\t*\n";
                out << "L\trep6\t+\trep7\t+\t*\n";
                out << "L\trep6\t+\trep8\t+\t*\n";
            }
        }

        // Block 3
        out << "S\t22\tGGGGG\n";
        out << "S\t23\tGGGGG\n";
        out << "S\t24\tGGGGG\n";
        out << "S\t25\tGGGGGTTTTT\n";
        // Before closure
        out << "L\t21\t+\t22\t+\t*\n";
        out << "L\t22\t+\t23\t+\t*\n";
        out << "L\t23\t+\t24\t+\t*\n";
        out << "L\t24\t+\t25\t+\t*\n";
        out << "L\t22\t+\t22\t+\t*\n";
        out << "L\t22\t+\t24\t+\t*\n";
        out << "L\t24\t+\t24\t+\t*\n";
        if (closure1) {
            // NOP
        }

        // Block 4
        out << "S\t26\tTTTTT\n";
        out << "S\t27\tTTTTT\n";
        out << "S\t28\tTTTTT\n";
        out << "S\t29\tTTTTTAAAAA\n";
        // Before closure
        out << "L\t25\t+\t26\t+\t*\n";
        out << "L\t26\t+\t27\t+\t*\n";
        out << "L\t27\t+\t28\t+\t*\n";
        out << "L\t28\t+\t29\t+\t*\n";
        out << "L\t25\t+\t26\t-\t*\n";
        out << "L\t26\t-\t27\t+\t*\n";
        out << "L\t26\t+\t28\t+\t*\n";
        out << "L\t27\t+\t28\t-\t*\n";
        out << "L\t28\t-\t29\t+\t*\n";
        if (closure1) {
            out << "L\t26\t+\t28\t-\t*\n";
            out << "L\t28\t-\t26\t+\t*\n";
            if (closure2) {
                out << "L\t26\t-\t28\t-\t*\n";
            }
        }

        // REP section
        // Block 1
        out << "S\t30\tAAAAA\n";
        out << "S\trep9\tCCCCC\n";
        out << "S\trep10\tGGGGG\n";
        out << "S\t31\tAAAAA\n";
        out << "S\t32\tAAAAA\n";
        out << "S\t33\tAAAAACCCCC\n";
        // Before closure
        out << "L\t29\t+\t30\t+\t*\n";
        out << "L\t30\t+\t31\t+\t*\n";
        out << "L\t31\t+\t32\t+\t*\n";
        out << "L\t32\t+\t33\t+\t*\n";
        out << "L\t29\t+\t30\t-\t*\n";
        out << "L\t30\t-\t31\t+\t*\n";
        out << "L\t30\t+\trep9\t+\t*\n";
        out << "L\t30\t+\trep10\t+\t*\n";
        out << "L\trep9\t+\t32\t+\t*\n";
        out << "L\trep10\t+\t32\t+\t*\n";
        out << "L\t31\t+\t32\t-\t*\n";
        out << "L\t32\t-\t33\t+\t*\n";
        if (closure1) {
            out << "L\t30\t-\trep9\t+\t*\n";
            out << "L\t30\t-\trep10\t+\t*\n";
            out << "L\t32\t+\trep9\t-\t*\n";
            out << "L\t32\t+\trep10\t-\t*\n";
            if (closure2) {
                // NOP
            }
        }

        // Block 2
        out << "S\t34\tCCCCC\n";
        out << "S\t35\tCCCCC\n";
        out << "S\trep11\tAAAAA\n";
        out << "S\trep12\tGGGGG\n";
        out << "S\t36\tCCCCC\n";
        out << "S\t37\tCCCCCGGGGG\n";
        // Before closure
        out << "L\t33\t+\t34\t+\t*\n";
        out << "L\t34\t+\t35\t+\t*\n";
        out << "L\t35\t+\t36\t+\t*\n";
        out << "L\t36\t+\t37\t+\t*\n";
        out << "L\t34\t+\t34\t+\t*\n";
        out << "L\t34\t+\trep11\t+\t*\n";
        out << "L\t34\t+\trep12\t+\t*\n";
        out << "L\trep11\t+\t36\t+\t*\n";
        out << "L\trep12\t+\t36\t+\t*\n";
        out << "L\t36\t+\t36\t+\t*\n";
        if (closure1) {
            // NOP
        }

        // DUP section
        // Block 1
        out << "S\t38\tGGGGG\n";
        out << "S\t39\tGGGGG\n";
        out << "S\t40\tGGGGG\n";
        out << "S\t41\tGGGGGTTTTT\n";
        // Before closure
        out << "L\t37\t+\t38\t+\t*\n";
        out << "L\t38\t+\t39\t+\t*\n";
        out << "L\t39\t+\t40\t+\t*\n";
        out << "L\t40\t+\t41\t+\t*\n";
        out << "L\t37\t+\t38\t-\t*\n";
        out << "L\t38\t-\t39\t+\t*\n";
        out << "L\t39\t+\t39\t+\t*\n";
        out << "L\t39\t+\t40\t-\t*\n";
        out << "L\t40\t-\t41\t+\t*\n";
        if (closure1) {
            // NOP
        }

        // Block 2
        out << "S\t42\tTTTTT\n";
        out << "S\t43\tTTTTT\n";
        out << "S\t44\tTTTTT\n";
        out << "S\t45\tTTTTTAAAAA\n";
        // Before closure
        out << "L\t41\t+\t42\t+\t*\n";
        out << "L\t42\t+\t43\t+\t*\n";
        out << "L\t43\t+\t44\t+\t*\n";
        out << "L\t44\t+\t45\t+\t*\n";
        out << "L\t42\t+\t42\t+\t*\n";
        out << "L\t43\t+\t43\t+\t*\n";
        out << "L\t44\t+\t44\t+\t*\n";
        if (closure1) {
            // NOP
        }

        // INV section
        // Block 1
        out << "S\t46\tAAAAA\n";
        out << "S\t47\tAAAAA\n";
        out << "S\t48\tAAAAA\n";
        out << "S\t49\tAAAAA\n";
        // Before closure
        out << "L\t45\t+\t46\t+\t*\n";
        out << "L\t46\t+\t47\t+\t*\n";
        out << "L\t47\t+\t48\t+\t*\n";
        out << "L\t48\t+\t49\t+\t*\n";
        out << "L\t45\t+\t46\t-\t*\n";
        out << "L\t46\t-\t47\t+\t*\n";
        out << "L\t46\t+\t47\t-\t*\n";
        out << "L\t47\t-\t48\t+\t*\n";
        out << "L\t47\t+\t48\t-\t*\n";
        out << "L\t48\t-\t49\t+\t*\n";
        if (closure1) {
            out << "L\t45\t+\t48\t+\t*\n";
            out << "L\t46\t-\t47\t-\t*\n";
            out << "L\t46\t+\t49\t+\t*\n";
            out << "L\t47\t-\t48\t-\t*\n";
            if (closure2) {
                out << "L\t45\t+\t48\t-\t*\n";
                out << "L\t46\t-\t49\t+\t*\n";
            }
        }
    }
    else if (dataset==1) {
        // REP section
        // Block 1
        out << "S\t1\tCCCCCCCCCCCCCCC\n";
        out << "S\t2\tCCCCC\n";
        out << "S\t3\tCCCCC\n";
        out << "S\t4\tCCCCC\n";
        out << "S\t5\tCCCCCCCCCCCCCCCGGGGG\n";
        out << "S\trep1\tAAAAA\n";
        out << "S\trep2\tGGGGG\n";
        out << "S\trep3\tAAAAA\n";
        out << "S\trep4\tGGGGG\n";
        out << "S\trep5\tAAAAA\n";
        out << "S\trep6\tGGGGG\n";
        // Before closure
        out << "L\t1\t+\t2\t+\t*\n";
        out << "L\t2\t+\t3\t+\t*\n";
        out << "L\t3\t+\t4\t+\t*\n";
        out << "L\t4\t+\t5\t+\t*\n";
        out << "L\t1\t+\trep1\t+\t*\n";
        out << "L\t1\t+\trep2\t+\t*\n";
        out << "L\trep1\t+\t3\t+\t*\n";
        out << "L\trep2\t+\t3\t+\t*\n";
        out << "L\t2\t+\trep3\t+\t*\n";
        out << "L\t2\t+\trep4\t+\t*\n";
        out << "L\trep3\t+\t4\t+\t*\n";
        out << "L\trep4\t+\t4\t+\t*\n";
        out << "L\t3\t+\trep5\t+\t*\n";
        out << "L\t3\t+\trep6\t+\t*\n";
        out << "L\trep5\t+\t5\t+\t*\n";
        out << "L\trep6\t+\t5\t+\t*\n";
        if (closure1) {
            out << "L\trep1\t+\trep3\t+\t*\n";
            out << "L\trep1\t+\trep4\t+\t*\n";
            out << "L\trep2\t+\trep3\t+\t*\n";
            out << "L\trep2\t+\trep4\t+\t*\n";
            out << "L\trep3\t+\trep5\t+\t*\n";
            out << "L\trep3\t+\trep6\t+\t*\n";
            out << "L\trep4\t+\trep5\t+\t*\n";
            out << "L\trep4\t+\trep6\t+\t*\n";
            if (closure2) {
                // NOP
            }
        }

        // DEL section
        // Block 1
        out << "S\t6\tGGGGG\n";
        out << "S\t7\tGGGGGTTTTT\n";
        out << "S\tins1\tCCCCC\n";
        out << "S\tins2\tTTTTT\n";
        out << "S\tins3\tCCCCC\n";
        out << "S\tins4\tTTTTT\n";
        // Before closure
        out << "L\t5\t+\t6\t+\t*\n";
        out << "L\t6\t+\t7\t+\t*\n";
        out << "L\t5\t+\tins1\t+\t*\n";
        out << "L\t5\t+\tins2\t+\t*\n";
        out << "L\tins1\t+\t6\t+\t*\n";
        out << "L\tins2\t+\t6\t+\t*\n";
        out << "L\t6\t+\tins3\t+\t*\n";
        out << "L\t6\t+\tins4\t+\t*\n";
        out << "L\tins3\t+\t7\t+\t*\n";
        out << "L\tins4\t+\t7\t+\t*\n";
        out << "L\t5\t+\t7\t+\t*\n";
        if (closure1) {
            out << "L\tins1\t+\t7\t+\t*\n";
            out << "L\tins2\t+\t7\t+\t*\n";
            out << "L\t5\t+\tins3\t+\t*\n";
            out << "L\t5\t+\tins4\t+\t*\n";
            if (closure2) {
                out << "L\tins1\t+\tins3\t+\t*\n";
                out << "L\tins1\t+\tins4\t+\t*\n";
                out << "L\tins2\t+\tins3\t+\t*\n";
                out << "L\tins2\t+\tins4\t+\t*\n";
            }
        }

        // REP section
        // Block 1
        out << "S\t8\tTTTTT\n";
        out << "S\t9\tTTTTTAAAAA\n";
        out << "S\tins5\tAAAAA\n";
        out << "S\tins6\tCCCCC\n";
        out << "S\trep7\tAAAAA\n";
        out << "S\trep8\tCCCCC\n";
        out << "S\tins7\tAAAAA\n";
        out << "S\tins8\tCCCCC\n";
        // Before closure
        out << "L\t7\t+\t8\t+\t*\n";
        out << "L\t8\t+\t9\t+\t*\n";
        out << "L\t7\t+\tins5\t+\t*\n";
        out << "L\t7\t+\tins6\t+\t*\n";
        out << "L\tins5\t+\t8\t+\t*\n";
        out << "L\tins6\t+\t8\t+\t*\n";
        out << "L\t7\t+\trep7\t+\t*\n";
        out << "L\t7\t+\trep8\t+\t*\n";
        out << "L\trep7\t+\t9\t+\t*\n";
        out << "L\trep8\t+\t9\t+\t*\n";
        out << "L\t8\t+\tins7\t+\t*\n";
        out << "L\t8\t+\tins8\t+\t*\n";
        out << "L\tins7\t+\t9\t+\t*\n";
        out << "L\tins8\t+\t9\t+\t*\n";
        // After closure
        if (closure1) {
            out << "L\tins5\t+\trep7\t+\t*\n";
            out << "L\tins5\t+\trep8\t+\t*\n";
            out << "L\tins6\t+\trep7\t+\t*\n";
            out << "L\tins6\t+\trep8\t+\t*\n";
            out << "L\trep7\t+\tins7\t+\t*\n";
            out << "L\trep7\t+\tins8\t+\t*\n";
            out << "L\trep8\t+\tins7\t+\t*\n";
            out << "L\trep8\t+\tins8\t+\t*\n";
            if (closure2) {
                // NOP
            }
        }

        // Block2
        out << "S\t10\tAAAAA\n";
        out << "S\t11\tAAAAA\n";
        out << "S\t12\tAAAAA\n";
        out << "S\t13\tAAAAACCCCC\n";
        out << "S\trep9\tCCCCC\n";
        out << "S\trep10\tGGGGG\n";
        // Before closure
        out << "L\t9\t+\t10\t+\t*\n";
        out << "L\t10\t+\t11\t+\t*\n";
        out << "L\t11\t+\t12\t+\t*\n";
        out << "L\t12\t+\t13\t+\t*\n";
        out << "L\t9\t+\t11\t+\t*\n";
        out << "L\t11\t+\t13\t+\t*\n";
        out << "L\t10\t+\trep9\t+\t*\n";
        out << "L\t10\t+\trep10\t+\t*\n";
        out << "L\trep9\t+\t12\t+\t*\n";
        out << "L\trep10\t+\t12\t+\t*\n";
        // After closure
        if (closure1) {
            out << "L\t9\t+\trep9\t+\t*\n";
            out << "L\t9\t+\trep10\t+\t*\n";
            out << "L\trep9\t+\t13\t+\t*\n";
            out << "L\trep10\t+\t13\t+\t*\n";
            if (closure2) {
                // NOP
            }
        }

        // INV section
        // Block 1
        out << "S\t14\tCCCCC\n";
        out << "S\t15\tCCCCCGGGGG\n";
        out << "S\tins9\tGGGGG\n";
        out << "S\tins10\tTTTTT\n";
        out << "S\tins11\tGGGGG\n";
        out << "S\tins12\tTTTTT\n";
        // Before closure
        out << "L\t13\t+\t14\t+\t*\n";
        out << "L\t14\t+\t15\t+\t*\n";
        out << "L\t13\t+\tins9\t+\t*\n";
        out << "L\t13\t+\tins10\t+\t*\n";
        out << "L\tins9\t+\t14\t+\t*\n";
        out << "L\tins10\t+\t14\t+\t*\n";
        out << "L\t13\t+\t14\t-\t*\n";
        out << "L\t14\t-\t15\t+\t*\n";
        out << "L\t14\t+\tins11\t+\t*\n";
        out << "L\t14\t+\tins12\t+\t*\n";
        out << "L\tins11\t+\t15\t+\t*\n";
        out << "L\tins12\t+\t15\t+\t*\n";
        // After closure
        if (closure1) {
            out << "L\t14\t+\tins9\t-\t*\n";
            out << "L\t14\t+\tins10\t-\t*\n";
            out << "L\t14\t-\tins11\t+\t*\n";
            out << "L\t14\t-\tins12\t+\t*\n";
            out << "L\t15\t-\tins9\t+\t*\n";
            out << "L\t15\t-\tins10\t+\t*\n";
            out << "L\tins11\t+\t13\t-\t*\n";
            out << "L\tins12\t+\t13\t-\t*\n";
            if (closure2) {
                out << "L\tins9\t+\tins11\t-\t*\n";
                out << "L\tins9\t+\tins12\t-\t*\n";
                out << "L\tins10\t+\tins11\t-\t*\n";
                out << "L\tins10\t+\tins12\t-\t*\n";
                out << "L\tins9\t-\tins11\t+\t*\n";
                out << "L\tins9\t-\tins12\t+\t*\n";
                out << "L\tins10\t-\tins11\t+\t*\n";
                out << "L\tins10\t-\tins12\t+\t*\n";
            }
        }

        // Block 2
        out << "S\t16\tGGGGG\n";
        out << "S\t17\tGGGGG\n";
        out << "S\t18\tGGGGG\n";
        out << "S\t19\tGGGGGTTTTT\n";
        // Before closure
        out << "L\t15\t+\t16\t+\t*\n";
        out << "L\t16\t+\t17\t+\t*\n";
        out << "L\t17\t+\t18\t+\t*\n";
        out << "L\t18\t+\t19\t+\t*\n";
        out << "L\t15\t+\t17\t+\t*\n";
        out << "L\t17\t+\t19\t+\t*\n";
        out << "L\t16\t+\t17\t-\t*\n";
        out << "L\t17\t-\t18\t+\t*\n";
        // After closure
        if (closure1) {
            out << "L\t15\t+\t17\t-\t*\n";
            out << "L\t17\t-\t19\t+\t*\n";
            if (closure2) {
                // NOP
            }
        }

        // Block 3
        out << "S\t20\tTTTTT\n";
        out << "S\t21\tTTTTT\n";
        out << "S\t22\tTTTTT\n";
        out << "S\t23\tTTTTTAAAAA\n";
        out << "S\trep11\tAAAAA\n";
        out << "S\trep12\tCCCCC\n";
        out << "S\trep13\tAAAAA\n";
        out << "S\trep14\tCCCCC\n";
        // Before closure
        out << "L\t19\t+\t20\t+\t*\n";
        out << "L\t20\t+\t21\t+\t*\n";
        out << "L\t21\t+\t22\t+\t*\n";
        out << "L\t22\t+\t23\t+\t*\n";
        out << "L\t19\t+\trep11\t+\t*\n";
        out << "L\t19\t+\trep12\t+\t*\n";
        out << "L\trep11\t+\t21\t+\t*\n";
        out << "L\trep12\t+\t21\t+\t*\n";
        out << "L\t21\t+\trep13\t+\t*\n";
        out << "L\t21\t+\trep14\t+\t*\n";
        out << "L\trep13\t+\t23\t+\t*\n";
        out << "L\trep14\t+\t23\t+\t*\n";
        out << "L\t20\t+\t21\t-\t*\n";
        out << "L\t21\t-\t22\t+\t*\n";
        // After closure
        if (closure1) {
            out << "L\trep11\t+\t21\t-\t*\n";
            out << "L\trep12\t+\t21\t-\t*\n";
            out << "L\t21\t-\trep13\t+\t*\n";
            out << "L\t21\t-\trep14\t+\t*\n";
            if (closure2) {
                // NOP
            }
        }

        // Block 4
        out << "S\t24\tAAAAA\n";
        out << "S\t25\tAAAAA\n";
        out << "S\t26\tAAAAA\n";
        out << "S\t27\tAAAAA\n";
        // Before closure
        out << "L\t23\t+\t24\t+\t*\n";
        out << "L\t24\t+\t25\t+\t*\n";
        out << "L\t25\t+\t26\t+\t*\n";
        out << "L\t26\t+\t27\t+\t*\n";
        out << "L\t24\t+\t24\t+\t*\n";
        out << "L\t24\t+\t25\t-\t*\n";
        out << "L\t25\t-\t26\t+\t*\n";
        out << "L\t26\t+\t26\t+\t*\n";
        // After closure
        if (closure1) {
            // NOP
            if (closure2) {
                // NOP
            }
        }
    }
}


/**
 * Remark: only feasible bidirected walks are tested, and only a subset of all possible ones.
 *
 * @param path_id >0; prints a selected path, whose node IDs come from the graph built by `VariantGraph`, not from the
 * graph built by `print_truth_gfa()`;
 * @param supported_records the procedure sets this to the IDs of the VCF records supported by `path_id`, in arbitrary
 * order.
 */
void load_true_path(size_t dataset, VariantGraph& graph, size_t path_id, const vector<string>& node_ids, vector<string>& supported_records, string& buffer) {
    string ins_id, ins_node, rep_id, rep_node, inv_id, inv_node;
    string left_rep_id, left_rep_node, right_rep_id, right_rep_node;
    string left_inv_id, left_inv_node, right_inv_id, right_inv_node;
    string left_ins_id, left_ins_node, right_ins_id, right_ins_node;

    graph.destroy_paths();
    supported_records.clear();

    if (dataset==0) {
        // INS section
        // Block 1
        ins_id="ins1"; ins_node="1";
        if (path_id==1) {
            graph.load_gfa_path("21+,"+ins_node+"+,23+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del1");
            supported_records.emplace_back(ins_id);
        }
        else if (path_id==2) {
            graph.load_gfa_path("22+,"+ins_node+"+,24+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(ins_id);
            supported_records.emplace_back("del2");
        }
        else if (path_id==3) {
            graph.load_gfa_path("21+,24+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del1");
            supported_records.emplace_back("del2");
        }
        else if (path_id==4) {
            graph.load_gfa_path("21+,"+ins_node+"+,24+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del1");
            supported_records.emplace_back(ins_id);
            supported_records.emplace_back("del2");
        }
        else if (path_id==5) {
            graph.load_gfa_path("21+,"+ins_node+"+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del1");
        }
        else if (path_id==6) {
            graph.load_gfa_path(ins_node+"+,24+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del2");
        }

        // Block 2
        left_rep_id="rep1"; left_rep_node="3";
        right_rep_id="rep3"; right_rep_node="7";
        ins_id="ins3"; ins_node="5";
        if (path_id==7) {
            graph.load_gfa_path("24+,"+left_rep_node+"+,26+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_rep_id);
        }
        else if (path_id==8) {
            graph.load_gfa_path("24+,"+left_rep_node+"+,"+ins_node+"+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_rep_id);
        }
        else if (path_id==9) {
            graph.load_gfa_path("24+,"+left_rep_node+"+,"+right_rep_node+"+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_rep_id);
        }
        else if (path_id==10) {
            graph.load_gfa_path("24+,"+left_rep_node+"+,"+ins_node+"+,26+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_rep_id);
            supported_records.emplace_back(ins_id);
        }
        else if (path_id==11) {
            graph.load_gfa_path("24+,"+left_rep_node+"+,"+ins_node+"+,"+right_rep_node+"+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_rep_id);
            supported_records.emplace_back(ins_id);
        }
        else if (path_id==12) {
            graph.load_gfa_path("24+,"+left_rep_node+"+,"+ins_node+"+,"+right_rep_node+"+,27+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_rep_id);
            supported_records.emplace_back(ins_id);
            supported_records.emplace_back(right_rep_id);
        }
        else if (path_id==13) {
            graph.load_gfa_path("24+,"+left_rep_node+"+,"+right_rep_node+"+,27+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_rep_id);
            supported_records.emplace_back(right_rep_id);
        }
        else if (path_id==14) {
            graph.load_gfa_path("25+,"+ins_node+"+,"+right_rep_node+"+,27+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(ins_id);
            supported_records.emplace_back(right_rep_id);
        }
        else if (path_id==15) {
            graph.load_gfa_path(left_rep_node+"+,"+ins_node+"+,"+right_rep_node+"+,27+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(ins_id);
            supported_records.emplace_back(right_rep_id);
        }
        else if (path_id==16) {
            graph.load_gfa_path(ins_node+"+,"+right_rep_node+"+,27+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(right_rep_id);
        }
        else if (path_id==17) {
            graph.load_gfa_path(left_rep_node+"+,"+right_rep_node+"+,27+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(right_rep_id);
        }
        else if (path_id==18) {
            graph.load_gfa_path("25+,"+right_rep_node+"+,27+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(right_rep_id);
        }

        // Block 3
        // NOP

        // Block 4
        ins_id="ins7"; ins_node="11"; left_inv_id="inv1"; left_inv_node="31"; right_inv_id="inv2"; right_inv_node="32";
        if (path_id==19) {
            graph.load_gfa_path("30+,"+left_inv_node+"-,32+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_inv_id);
        }
        else if (path_id==20) {
            graph.load_gfa_path("30+,"+left_inv_node+"-,"+ins_node+"+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_inv_id);
        }
        else if (path_id==21) {
            graph.load_gfa_path("30+,"+left_inv_node+"-,32-",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_inv_id);
        }
        else if (path_id==22) {
            graph.load_gfa_path("30+,"+left_inv_node+"-,"+ins_node+"+,32+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_inv_id);
            supported_records.emplace_back(ins_id);
        }
        else if (path_id==23) {
            graph.load_gfa_path("30+,"+left_inv_node+"-,"+ins_node+"+,32-",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_inv_id);
            supported_records.emplace_back(ins_id);
        }
        else if (path_id==24) {
            graph.load_gfa_path("30+,"+left_inv_node+"-,"+ins_node+"+,30-",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_inv_id);
            supported_records.emplace_back(ins_id);
        }
        else if (path_id==25) {
            graph.load_gfa_path("30+,"+left_inv_node+"-,"+ins_node+"+,"+right_inv_node+"-,33+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_inv_id);
            supported_records.emplace_back(ins_id);
            supported_records.emplace_back(right_inv_id);
        }
        else if (path_id==26) {
            graph.load_gfa_path("30+,"+left_inv_node+"-,"+right_inv_node+"-,33+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_inv_id);
            supported_records.emplace_back(right_inv_id);
        }
        else if (path_id==27) {
            graph.load_gfa_path("31+,"+ins_node+"+,32+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(ins_id);
        }
        else if (path_id==28) {
            graph.load_gfa_path(left_inv_node+"+,"+ins_node+"+,"+right_inv_node+"-,33+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(ins_id);
            supported_records.emplace_back(right_inv_id);
        }
        else if (path_id==29) {
            graph.load_gfa_path(left_inv_node+"-,"+ins_node+"+,"+right_inv_node+"-,33+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(ins_id);
            supported_records.emplace_back(right_inv_id);
        }
        else if (path_id==30) {
            graph.load_gfa_path("33-,"+ins_node+"+,"+right_inv_node+"-,33+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(ins_id);
            supported_records.emplace_back(right_inv_id);
        }
        else if (path_id==31) {
            graph.load_gfa_path(ins_node+"+,"+right_inv_node+"-,33+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(right_inv_id);
        }
        else if (path_id==32) {
            graph.load_gfa_path(left_inv_node+"+,"+right_inv_node+"-,33+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(right_inv_id);
        }
        else if (path_id==33) {
            graph.load_gfa_path(left_inv_node+"-,"+right_inv_node+"-,33+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(right_inv_id);
        }
        else if (path_id==34) {
            graph.load_gfa_path("30+,"+ins_node+"-,"+left_inv_node+"+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(ins_id);
            supported_records.emplace_back(left_inv_id);
        }
        else if (path_id==35) {
            graph.load_gfa_path(right_inv_node+"+,"+ins_node+"-,"+left_inv_node+"+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(ins_id);
        }
        else if (path_id==36) {
            graph.load_gfa_path("30+,"+ins_node+"-,33+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(ins_id);
        }
        else if (path_id==37) {
            graph.load_gfa_path(right_inv_node+"+,"+ins_node+"-,33+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(right_inv_id);
            supported_records.emplace_back(ins_id);
        }
        else if (path_id==38) {
            graph.load_gfa_path(left_inv_node+"-,"+right_inv_node+"-",node_ids,to_string(path_id),buffer);
        }
        else if (path_id==39) {
            graph.load_gfa_path("30+,33+",node_ids,to_string(path_id),buffer);
        }

        // DEL section
        // Block 1
        if (path_id==40) {
            graph.load_gfa_path("33+,36+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del3");
            supported_records.emplace_back("del4");
        }
        else if (path_id==41) {
            graph.load_gfa_path("34+,37+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del4");
            supported_records.emplace_back("del5");
        }
        else if (path_id==42) {
            graph.load_gfa_path("33+,37+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del3");
            supported_records.emplace_back("del4");
            supported_records.emplace_back("del5");
        }

        // Block 2
        left_rep_id="rep5"; left_rep_node="13"; right_rep_id="rep7"; right_rep_node="15";
        if (path_id==43) {
            graph.load_gfa_path("37+,"+left_rep_node+"+,39+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_rep_id);
        }
        else if (path_id==44) {
            graph.load_gfa_path("37+,"+left_rep_node+"+,40+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_rep_id);
            supported_records.emplace_back("del6");
        }
        else if (path_id==45) {
            graph.load_gfa_path("37+,"+left_rep_node+"+,"+right_rep_node+"+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_rep_id);
            supported_records.emplace_back("del6");
        }
        else if (path_id==46) {
            graph.load_gfa_path("37+,"+left_rep_node+"+,"+right_rep_node+"+,41+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_rep_id);
            supported_records.emplace_back("del6");
            supported_records.emplace_back(right_rep_id);
        }
        else if (path_id==47) {
            graph.load_gfa_path("38+,40+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del6");
        }
        else if (path_id==48) {
            graph.load_gfa_path("38+,"+right_rep_node+"+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del6");
        }
        else if (path_id==49) {
            graph.load_gfa_path(left_rep_node+"+,40+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del6");
        }
        else if (path_id==50) {
            graph.load_gfa_path(left_rep_node+"+,"+right_rep_node+"+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del6");
        }
        else if (path_id==51) {
            graph.load_gfa_path("38+,"+right_rep_node+"+,41+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del6");
            supported_records.emplace_back(right_rep_id);
        }
        else if (path_id==52) {
            graph.load_gfa_path(left_rep_node+"+,"+right_rep_node+"+,41+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del6");
            supported_records.emplace_back(right_rep_id);
        }
        else if (path_id==53) {
            graph.load_gfa_path("39+,"+right_rep_node+"+,41+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(right_rep_id);
        }

        // Block 3
        // NOP

        // Block 4
        left_inv_id="inv3"; left_inv_node="46"; right_inv_id="inv4"; right_inv_node="48";
        if (path_id==54) {
            graph.load_gfa_path("45+,"+left_inv_node+"-,47+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_inv_id);
        }
        else if (path_id==55) {
            graph.load_gfa_path("45+,"+left_inv_node+"-,"+right_inv_node+"+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_inv_id);
            supported_records.emplace_back("del8");
        }
        else if (path_id==56) {
            graph.load_gfa_path("45+,"+left_inv_node+"-,"+right_inv_node+"-,49+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_inv_id);
            supported_records.emplace_back("del8");
            supported_records.emplace_back(right_inv_id);
        }
        else if (path_id==57) {
            graph.load_gfa_path(left_inv_node+"+,"+right_inv_node+"+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del8");
        }
        else if (path_id==58) {
            graph.load_gfa_path(right_inv_node+"+,"+left_inv_node+"+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del8");
        }
        else if (path_id==59) {
            graph.load_gfa_path(left_inv_node+"+,"+right_inv_node+"-,"+left_inv_node+"+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del8");
        }
        else if (path_id==60) {
            graph.load_gfa_path(left_inv_node+"+,"+right_inv_node+"-,49+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del8");
            supported_records.emplace_back(right_inv_id);
        }
        else if (path_id==61) {
            graph.load_gfa_path("47+,"+right_inv_node+"-,49+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(right_inv_id);
        }

        // REP section
        // Block 1
        left_inv_id="inv5"; left_inv_node="50"; right_inv_id="inv6"; right_inv_node="52"; rep_id="rep9"; rep_node="17";
        if (path_id==62) {
            graph.load_gfa_path("49+,"+left_inv_node+"-,"+rep_node+"+,"+right_inv_node+"+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_inv_id);
            supported_records.emplace_back(rep_id);
        }
        else if (path_id==63) {
            graph.load_gfa_path("49+,"+left_inv_node+"-,"+rep_node+"+,"+right_inv_node+"-,53+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_inv_id);
            supported_records.emplace_back(rep_id);
            supported_records.emplace_back(right_inv_id);
        }
        else if (path_id==64) {
            graph.load_gfa_path(left_inv_node+"+,"+rep_node+"+,"+right_inv_node+"+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(rep_id);
        }
        else if (path_id==65) {
            graph.load_gfa_path(left_inv_node+"-,"+rep_node+"+,"+right_inv_node+"-",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(rep_id);
        }
        else if (path_id==66) {
            graph.load_gfa_path(left_inv_node+"+,"+rep_node+"+,"+right_inv_node+"-,53+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(rep_id);
            supported_records.emplace_back(right_inv_id);
        }

        // Block 2
        // NOP

        // DUP section
        // Block 1
        // NOP

        // Block 2
        // NOP

        // INV section
        // Block 1
        left_inv_id="inv9"; left_inv_node="66"; inv_id="inv10"; inv_node="67"; right_inv_id="inv11"; right_inv_node="68";
        if (path_id==67) {
            graph.load_gfa_path("65+,"+left_inv_node+"-,"+inv_node+"-",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_inv_id);
        }
        else if (path_id==68) {
            graph.load_gfa_path("65+,"+left_inv_node+"-,69+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_inv_id);
        }
        else if (path_id==69) {
            graph.load_gfa_path("65+,"+left_inv_node+"-,"+inv_node+"-,"+right_inv_node+"+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_inv_id);
            supported_records.emplace_back(inv_id);
        }
        else if (path_id==70) {
            graph.load_gfa_path("65+,"+left_inv_node+"-,"+inv_node+"-,"+right_inv_node+"-,69+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_inv_id);
            supported_records.emplace_back(inv_id);
            supported_records.emplace_back(right_inv_id);
        }
        else if (path_id==71) {
            graph.load_gfa_path("65+,"+right_inv_node+"+,"+inv_node+"+,"+left_inv_node+"+,69+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_inv_id);
            supported_records.emplace_back(inv_id);
            supported_records.emplace_back(right_inv_id);
        }
        else if (path_id==72) {
            graph.load_gfa_path("65+,"+right_inv_node+"-,"+inv_node+"-,"+left_inv_node+"-,69+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_inv_id);
            supported_records.emplace_back(inv_id);
            supported_records.emplace_back(right_inv_id);
        }
        else if (path_id==73) {
            graph.load_gfa_path(left_inv_node+"+,"+inv_node+"-,"+right_inv_node+"-,69+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(inv_id);
            supported_records.emplace_back(right_inv_id);
        }
        else if (path_id==74) {
            graph.load_gfa_path(inv_node+"-,"+right_inv_node+"-,69+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(right_inv_id);
        }
        else if (path_id==75) {
            graph.load_gfa_path("65+,"+right_inv_node+"-,69+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(right_inv_id);
        }
        else if (path_id==76) {
            graph.load_gfa_path("65+,"+right_inv_node+"+",node_ids,to_string(path_id),buffer);
        }
        else if (path_id==77) {
            graph.load_gfa_path(left_inv_node+"+,69+",node_ids,to_string(path_id),buffer);
        }
        else if (path_id==78) {
            graph.load_gfa_path(left_inv_node+"+,"+inv_node+"+,"+left_inv_node+"+",node_ids,to_string(path_id),buffer);
        }
        else if (path_id==79) {
            graph.load_gfa_path(inv_node+"+,"+right_inv_node+"+,"+inv_node+"+",node_ids,to_string(path_id),buffer);
        }
    }
    else if (dataset==1) {
        // REP section
        // Block 1
        left_rep_id="rep1"; left_rep_node="1"; rep_id="rep3"; rep_node="3"; right_rep_id="rep5"; right_rep_node="5";
        if (path_id==1) {
            graph.load_gfa_path("27+,"+left_rep_node+"+,"+rep_node+"+,30+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_rep_id);
            supported_records.emplace_back(rep_id);
        }
        else if (path_id==2) {
            graph.load_gfa_path("28+,"+rep_node+"+,"+right_rep_node+"+,31+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(rep_id);
            supported_records.emplace_back(right_rep_id);
        }
        else if (path_id==3) {
            graph.load_gfa_path("27+,"+left_rep_node+"+,"+rep_node+"+,"+right_rep_node+"+,31+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_rep_id);
            supported_records.emplace_back(rep_id);
            supported_records.emplace_back(right_rep_id);
        }
        else if (path_id==4) {
            graph.load_gfa_path("27+,"+left_rep_node+"+,"+rep_node+"+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_rep_id);
        }
        else if (path_id==5) {
            graph.load_gfa_path(left_rep_node+"+,"+rep_node+"+,"+right_rep_node+"+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(rep_id);
        }
        else if (path_id==6) {
            graph.load_gfa_path(rep_node+"+,"+right_rep_node+"+,31+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(right_rep_id);
        }

        // DEL section
        // Block 1
        left_ins_id="ins1"; left_ins_node="7"; right_ins_id="ins3"; right_ins_node="9";
        if (path_id==7) {
            graph.load_gfa_path("31+,"+left_ins_node+"+,33+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_ins_id);
            supported_records.emplace_back("del1");
        }
        else if (path_id==8) {
            graph.load_gfa_path("31+,"+left_ins_node+"+,"+right_ins_node+"+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_ins_id);
            supported_records.emplace_back("del1");
        }
        else if (path_id==9) {
            graph.load_gfa_path("31+,"+right_ins_node+"+,33+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del1");
            supported_records.emplace_back(right_ins_id);
        }
        else if (path_id==10) {
            graph.load_gfa_path(left_ins_node+"+,"+right_ins_node+"+,33+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del1");
            supported_records.emplace_back(right_ins_id);
        }
        else if (path_id==11) {
            graph.load_gfa_path("31+,"+left_ins_node+"+,"+right_ins_node+"+,33+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_ins_id);
            supported_records.emplace_back("del1");
            supported_records.emplace_back(right_ins_id);
        }
        else if (path_id==12) {
            graph.load_gfa_path(left_ins_node+"+,"+right_ins_node+"+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del1");
        }

        // REP section
        // Block 1
        left_ins_id="ins5"; left_ins_node="11"; rep_id="rep7"; rep_node="13"; right_ins_id="ins7"; right_ins_node="15";
        if (path_id==13) {
            graph.load_gfa_path("33+,"+left_ins_node+"+,"+rep_node+"+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_ins_id);
        }
        else if (path_id==14) {
            graph.load_gfa_path("33+,"+left_ins_node+"+,"+rep_node+"+,35+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_ins_id);
            supported_records.emplace_back(rep_id);
        }
        else if (path_id==15) {
            graph.load_gfa_path("33+,"+left_ins_node+"+,"+rep_node+"+,"+right_ins_node+"+,35+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_ins_id);
            supported_records.emplace_back(rep_id);
            supported_records.emplace_back(right_ins_id);
        }
        else if (path_id==16) {
            graph.load_gfa_path(left_ins_node+"+,"+rep_node+"+,"+right_ins_node+"+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(rep_id);
        }
        else if (path_id==17) {
            graph.load_gfa_path("33+,"+rep_node+"+,"+right_ins_node+"+,35+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(rep_id);
            supported_records.emplace_back(right_ins_id);
        }
        else if (path_id==18) {
            graph.load_gfa_path(rep_node+"+,"+right_ins_node+"+,35+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(right_ins_id);
        }

        // Block 2
        rep_id="rep9"; rep_node="17";
        if (path_id==19) {
            graph.load_gfa_path("35+,"+rep_node+"+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del2");
        }
        else if (path_id==20) {
            graph.load_gfa_path("35+,"+rep_node+"+,38+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del2");
            supported_records.emplace_back(rep_id);
        }
        else if (path_id==21) {
            graph.load_gfa_path("35+,"+rep_node+"+,39+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del2");
            supported_records.emplace_back(rep_id);
            supported_records.emplace_back("del3");
        }
        else if (path_id==22) {
            graph.load_gfa_path("36+,"+rep_node+"+,39+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(rep_id);
            supported_records.emplace_back("del3");
        }
        else if (path_id==23) {
            graph.load_gfa_path(rep_node+"+,39+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del3");
        }

        // INV section
        // Block 1
        left_ins_id="ins9"; left_ins_node="19"; right_ins_id="ins11"; right_ins_node="21";
        if (path_id==24) {
            graph.load_gfa_path("39+,"+left_ins_node+"+,40-",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_ins_id);
        }
        else if (path_id==25) {
            graph.load_gfa_path("39+,"+left_ins_node+"+,"+right_ins_node+"-",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_ins_id);
        }
        else if (path_id==26) {
            graph.load_gfa_path(right_ins_node+"-,"+left_ins_node+"+,40+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_ins_id);
        }
        else if (path_id==27) {
            graph.load_gfa_path("41-,"+left_ins_node+"+,"+right_ins_node+"-",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_ins_id);
            supported_records.emplace_back("inv12");
        }
        else if (path_id==28) {
            graph.load_gfa_path("39+,"+left_ins_node+"+,40-,41+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_ins_id);
            supported_records.emplace_back("inv12");
        }
        else if (path_id==29) {
            graph.load_gfa_path("40+,"+left_ins_node+"-,41+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_ins_id);
            supported_records.emplace_back("inv12");
        }
        else if (path_id==30) {
            graph.load_gfa_path(right_ins_node+"-,"+left_ins_node+"+,40-",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_ins_id);
            supported_records.emplace_back("inv12");
        }
        else if (path_id==31) {
            graph.load_gfa_path("41-,"+left_ins_node+"+,40-",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_ins_id);
            supported_records.emplace_back("inv12");
        }
        else if (path_id==32) {
            graph.load_gfa_path("39+,"+left_ins_node+"+,40-,"+right_ins_node+"+,41+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_ins_id);
            supported_records.emplace_back("inv12");
            supported_records.emplace_back(right_ins_id);
        }
        else if (path_id==33) {
            graph.load_gfa_path(left_ins_node+"+,"+right_ins_node+"-,"+left_ins_node+"+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_ins_id);
            supported_records.emplace_back("inv12");
            supported_records.emplace_back(right_ins_id);
        }
        else if (path_id==34) {
            graph.load_gfa_path(left_ins_node+"+,"+right_ins_node+"-,40+,"+left_ins_node+"-,"+right_ins_node+"+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_ins_id);
            supported_records.emplace_back("inv12");
            supported_records.emplace_back(right_ins_id);
        }
        else if (path_id==35) {
            graph.load_gfa_path("39+,"+right_ins_node+"-,40+,"+left_ins_node+"-,41+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_ins_id);
            supported_records.emplace_back("inv12");
            supported_records.emplace_back(right_ins_id);
        }
        else if (path_id==36) {
            graph.load_gfa_path("39+,40-,"+right_ins_node+"+,41+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("inv12");
            supported_records.emplace_back(right_ins_id);
        }
        else if (path_id==37) {
            graph.load_gfa_path("40-,"+right_ins_node+"+,39-",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("inv12");
            supported_records.emplace_back(right_ins_id);
        }
        else if (path_id==38) {
            graph.load_gfa_path("40-,"+right_ins_node+"+,"+left_ins_node+"-",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("inv12");
            supported_records.emplace_back(right_ins_id);
        }
        else if (path_id==39) {
            graph.load_gfa_path("39+,"+right_ins_node+"-,40+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("inv12");
            supported_records.emplace_back(right_ins_id);
        }
        else if (path_id==40) {
            graph.load_gfa_path("39+,"+right_ins_node+"-,"+left_ins_node+"+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("inv12");
            supported_records.emplace_back(right_ins_id);
        }
        else if (path_id==41) {
            graph.load_gfa_path("40-,"+right_ins_node+"+,41+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(right_ins_id);
        }
        else if (path_id==42) {
            graph.load_gfa_path(left_ins_node+"-,"+right_ins_node+"+,41+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(right_ins_id);
        }
        else if (path_id==43) {
            graph.load_gfa_path("40+,"+right_ins_node+"+,"+left_ins_node+"-",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(right_ins_id);
        }

        // Block 2
        else if (path_id==44) {
            graph.load_gfa_path("41+,43-",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del4");
        }
        else if (path_id==45) {
            graph.load_gfa_path("41+,43-,44+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del4");
            supported_records.emplace_back("inv13");
        }
        else if (path_id==46) {
            graph.load_gfa_path("41+,43-,45+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del4");
            supported_records.emplace_back("inv13");
            supported_records.emplace_back("del5");
        }
        else if (path_id==47) {
            graph.load_gfa_path("42+,43-,45+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("inv13");
            supported_records.emplace_back("del5");
        }
        else if (path_id==48) {
            graph.load_gfa_path("43-,45+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("del5");
        }

        // Block 3
        left_rep_id="rep11"; left_rep_node="23"; right_rep_id="rep13"; right_rep_node="25";
        if (path_id==49) {
            graph.load_gfa_path("45+,"+left_rep_node+"+,47-",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_rep_id);
        }
        else if (path_id==50) {
            graph.load_gfa_path("45+,"+left_rep_node+"+,47-,48+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_rep_id);
            supported_records.emplace_back("inv14");
        }
        else if (path_id==51) {
            graph.load_gfa_path("45+,"+left_rep_node+"+,47-,"+right_rep_node+"+,49+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(left_rep_id);
            supported_records.emplace_back("inv14");
            supported_records.emplace_back(right_rep_id);
        }
        else if (path_id==52) {
            graph.load_gfa_path(left_rep_node+"+,47-,"+right_rep_node+"+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("inv14");
        }
        else if (path_id==53) {
            graph.load_gfa_path("46+,47-,"+right_rep_node+"+,49+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back("inv14");
            supported_records.emplace_back(right_rep_id);
        }
        else if (path_id==54) {
            graph.load_gfa_path("47-,"+right_rep_node+"+,49+",node_ids,to_string(path_id),buffer);
            supported_records.emplace_back(right_rep_id);
        }

        // Block 4
        // NOP
    }
}


/**
 * @return X such that path IDs are in [1..X].
 */
size_t get_n_true_paths(size_t dataset) {
    if (dataset==0) return 79;
    else if (dataset==1) return 54;
    else return 0;
}


void print_gfa_colors(size_t dataset, ofstream& out) {
    const string COLOR_REF = "gray";
    const string COLOR_INS = "red";
    const string COLOR_REP = "green";

    int32_t i;

    out << "Name,Colour\n";
    if (dataset==0) {
        for (i=1; i<=49; i++) out << to_string(i) << "," << COLOR_REF << "\n";
        for (i=1; i<=8; i++) out << "ins" << to_string(i) << "," << COLOR_INS << "\n";
        for (i=1; i<=12; i++) out << "rep" << to_string(i) << "," << COLOR_REP << "\n";
    }
    else if (dataset==1) {
        for (i=1; i<=27; i++) out << to_string(i) << "," << COLOR_REF << "\n";
        for (i=1; i<=12; i++) out << "ins" << to_string(i) << "," << COLOR_INS << "\n";
        for (i=1; i<=14; i++) out << "rep" << to_string(i) << "," << COLOR_REP << "\n";
    }
}


unordered_map<string,string> get_chromosomes(size_t dataset) {
    unordered_map<string,string> out;
    if (dataset==0) out.emplace("chr1","AAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAA");
    else if (dataset==1) out.emplace("chr1","CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAA");
    return out;
}


/**
 * Zero-based, non-overlapping, [x..y).
 */
unordered_map<string,vector<interval_t>> get_tandem_track() {
    unordered_map<string,vector<interval_t>> out;
    return out;
}


int main(int argc, char* argv[]) {
    const path ROOT_DIR = path(argv[1]);
    const int32_t CLOSURE_1 = stoi(argv[2]);
    const int32_t CLOSURE_2 = stoi(argv[3]);

    const path INPUT_VCF = ROOT_DIR/"input.vcf";
    const path TRUTH_GFA = ROOT_DIR/"truth.gfa";
    const path TRUTH_GFA_COLORS = ROOT_DIR/"truth.csv";
    const path TEST_GFA = ROOT_DIR/"test.gfa";
    const int32_t SIGNATURE_N_STEPS = 10;
    const int32_t FLANK_LENGTH = 1000;
    const int32_t INTERIOR_FLANK_LENGTH = 1000;

    const unordered_map<string,vector<interval_t>> tandem_track = get_tandem_track();
    size_t i, j;
    size_t n_records, dataset, n_true_paths;
    string command, buffer;
    vector<string> supported_records, node_ids;
    vector<VcfRecord> records;
    unordered_map<string,string> chromosomes;

    for (dataset=0; dataset<=1; dataset++) {
        cerr << "Testing dataset " << to_string(dataset) << "...\n";

        // Printing truth files
        chromosomes=get_chromosomes(dataset);
        ofstream input_vcf(INPUT_VCF.string());
        print_truth_vcf_header(input_vcf);
        print_truth_vcf(dataset,input_vcf);
        input_vcf.close();
        ofstream truth_gfa(TRUTH_GFA.string());
        print_truth_gfa(dataset,truth_gfa,CLOSURE_1==1,CLOSURE_2==1);
        truth_gfa.close();
        ofstream truth_gfa_colors(TRUTH_GFA_COLORS.string());
        print_gfa_colors(dataset,truth_gfa_colors);
        truth_gfa_colors.close();

        // Loading the VCF and building the graph
        VcfReader reader(INPUT_VCF);
        reader.for_record_in_vcf([&](VcfRecord& record) {
            if ( (record.sv_type==VcfReader::TYPE_INSERTION && record.is_symbolic) ||
                 ((record.sv_type==VcfReader::TYPE_DELETION || record.sv_type==VcfReader::TYPE_INVERSION || record.sv_type==VcfReader::TYPE_DUPLICATION || record.sv_type==VcfReader::TYPE_REPLACEMENT) && record.sv_length==INT32_MAX)
               ) return;
            records.push_back(record);
        });
        n_records=records.size();
        VariantGraph graph(chromosomes,tandem_track);
        graph.build(records,FLANK_LENGTH,INTERIOR_FLANK_LENGTH);
        graph.to_gfa(TEST_GFA);

        cerr << "Testing to_gfa(): node sequences...\n";
        command.clear();
        command.append("grep ^S " + TRUTH_GFA.string() + " | cut -f 1,3 | sort > tmp1.txt");
        run_command(command);
        command.clear();
        command.append("grep ^S " + TEST_GFA.string() + " | cut -f 1,3 | sort > tmp2.txt");
        run_command(command);
        command.clear();
        command.append("diff --brief tmp1.txt tmp2.txt");
        run_command(command);

        cerr << "Testing to_gfa(): n. edges...\n";
        command.clear();
        command.append("grep ^L " + TRUTH_GFA.string() + " | wc -l > tmp1.txt");
        run_command(command);
        command.clear();
        command.append("grep ^L " + TEST_GFA.string() + " | wc -l > tmp2.txt");
        run_command(command);
        command.clear();
        command.append("diff --brief tmp1.txt tmp2.txt");
        run_command(command);

        cerr << "Testing to_gfa(): local topology (" << SIGNATURE_N_STEPS << " steps)...\n";
        graph.print_graph_signature(SIGNATURE_N_STEPS, "tmp1.txt");
        graph.load_gfa(TRUTH_GFA);
        graph.print_graph_signature(SIGNATURE_N_STEPS, "tmp2.txt");
        command.clear(); command.append("diff --brief tmp1.txt tmp2.txt"); run_command(command);

        cerr << "Testing edge-record correspondence...\n";
        n_true_paths=get_n_true_paths(dataset);
        for (i=1; i<=n_true_paths; i++) {
            cerr << "Testing path " << i << ":\n";
            load_true_path(dataset,graph,i,node_ids,supported_records,buffer);
            ofstream supported_truth_vcf("supported_truth.vcf");
            for (auto& id: supported_records) {
                for (j=0; j<n_records; j++) {
                    if (graph.vcf_records.at(j).id==id) {
                        graph.vcf_records.at(j).print(supported_truth_vcf);
                        supported_truth_vcf << '\n';
                        break;
                    }
                }
            }
            supported_truth_vcf.close();
            ofstream supported_test_vcf("supported_test.vcf"); ofstream unsupported_test_vcf("unsupported_test.vcf");
            graph.print_supported_vcf_records(supported_test_vcf,unsupported_test_vcf,false);
            supported_test_vcf.close(); unsupported_test_vcf.close();
            command.clear(); command.append("sort supported_test.vcf > supported_test_sorted.vcf"); run_command(command);
            command.clear(); command.append("sort supported_truth.vcf > supported_truth_sorted.vcf"); run_command(command);
            command.clear(); command.append("diff --brief supported_test_sorted.vcf supported_truth_sorted.vcf"); run_command(command);
        }

        cerr << "Removing temporary files...\n";
        command.clear(); command.append("rm -f tmp*.txt supported_*.vcf unsupported_*.vcf"); run_command(command);
    }
}