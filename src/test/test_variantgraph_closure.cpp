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


void print_truth_vcf(ofstream& out) {
    const string QUAL = "60";  // Arbitrary
    const string FILTER = "PASS";
    const string FORMAT = "GT";
    const string GT = "0/1";  // Arbitrary
    const string INFIX = QUAL+"\t"+FILTER;
    const string SUFFIX = FORMAT+"\t"+GT+"\n";

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


/**
 * @param path_id >0; prints a selected path;
 * @param supported_records the procedure sets this to the IDs of the VCF records supported by `path_id`, in arbitrary
 * order.
 */
void print_truth_gfa(ofstream& out, bool closure1, bool closure2, int32_t path_id, vector<string>& supported_records) {
    supported_records.clear();

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
    if (path_id==1) {
        out << "P\t" << to_string(path_id) << "\t1+,ins1+,3+\t*\n";
        supported_records.emplace_back("del1");
        supported_records.emplace_back("ins1");
    }
    if (path_id==2) {
        out << "P\t" << to_string(path_id) << "\t1+,ins1+,4+\t*\n";
        supported_records.emplace_back("del1");
        supported_records.emplace_back("ins1");
        supported_records.emplace_back("del2");
    }
    if (path_id==3) {
        out << "P\t" << to_string(path_id) << "\t1+,ins2+,3+\t*\n";
        supported_records.emplace_back("del1");
        supported_records.emplace_back("ins2");
    }
    if (path_id==4) {
        out << "P\t" << to_string(path_id) << "\t1+,ins2+,4+\t*\n";
        supported_records.emplace_back("del1");
        supported_records.emplace_back("ins2");
        supported_records.emplace_back("del2");
    }
    if (path_id==5) {
        out << "P\t" << to_string(path_id) << "\t1+,4+\t*\n";
        supported_records.emplace_back("del1");
        supported_records.emplace_back("del2");
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
    string left_rep_id = "rep1";
    string right_rep_id = "rep2";
    string ins_id = "ins3";  // ins3 or ins4
    if (path_id==6) {
        out << "P\t" << to_string(path_id) << "\t4+," << left_rep_id << "+,6+\t*\n";
        supported_records.emplace_back(left_rep_id);
    }
    if (path_id==7) {
        out << "P\t" << to_string(path_id) << "\t4+," << left_rep_id << "+," << ins_id << "+\t*\n";
        supported_records.emplace_back(left_rep_id);
        supported_records.emplace_back(ins_id);
    }
    if (path_id==8) {
        out << "P\t" << to_string(path_id) << "\t4+," << left_rep_id << "+," << ins_id << "+,6+\t*\n";
        supported_records.emplace_back(left_rep_id);
        supported_records.emplace_back(ins_id);
    }
    if (path_id==9) {
        out << "P\t" << to_string(path_id) << "\t4+," << left_rep_id << "+," << ins_id << "+," << right_rep_id << "+\t*\n";
        supported_records.emplace_back(left_rep_id);
        supported_records.emplace_back(ins_id);
        supported_records.emplace_back(right_rep_id);
    }
    if (path_id==10) {
        out << "P\t" << to_string(path_id) << "\t4+," << left_rep_id << "+," << ins_id << "+," << right_rep_id << "+,7+\t*\n";
        supported_records.emplace_back(left_rep_id);
        supported_records.emplace_back(ins_id);
        supported_records.emplace_back(right_rep_id);
    }
    if (path_id==11) {
        out << "P\t" << to_string(path_id) << "\t4+," << left_rep_id << "+," << right_rep_id << "+\t*\n";
        supported_records.emplace_back(left_rep_id);
        supported_records.emplace_back(right_rep_id);
    }
    if (path_id==12) {
        out << "P\t" << to_string(path_id) << "\t4+," << left_rep_id << "+," << right_rep_id << "+,7+\t*\n";
        supported_records.emplace_back(left_rep_id);
        supported_records.emplace_back(right_rep_id);
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
    ins_id="ins7";
    if (path_id==13) {
        out << "P\t" << to_string(path_id) << "\t10+,11-,12-,13+\t*\n";
        supported_records.emplace_back("inv1");
        supported_records.emplace_back("inv2");
    }
    if (path_id==14) {
        out << "P\t" << to_string(path_id) << "\t10+,13+\t*\n";
        supported_records.emplace_back("inv1");
        supported_records.emplace_back("inv2");
    }
    if (path_id==15) {
        out << "P\t" << to_string(path_id) << "\t10+,11-," << ins_id << "+,12-,13+\t*\n";
        supported_records.emplace_back("inv1");
        supported_records.emplace_back("inv2");
        supported_records.emplace_back(ins_id);
    }
    if (path_id==16) {
        out << "P\t" << to_string(path_id) << "\t10+," << ins_id << "-,13+\t*\n";
        supported_records.emplace_back("inv1");
        supported_records.emplace_back("inv2");
        supported_records.emplace_back(ins_id);
    }
    if (path_id==17) {
        out << "P\t" << to_string(path_id) << "\t10+,11-," << ins_id << "+,12+,13+\t*\n";
        supported_records.emplace_back("inv1");
        supported_records.emplace_back(ins_id);
    }
    if (path_id==18) {
        out << "P\t" << to_string(path_id) << "\t10+,11+," << ins_id << "+,12-,13+\t*\n";
        supported_records.emplace_back("inv2");
        supported_records.emplace_back(ins_id);
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
    if (path_id==19) {
        out << "P\t" << to_string(path_id) << "\t13+,16+\t*\n";
        supported_records.emplace_back("del3");
        supported_records.emplace_back("del4");
    }
    if (path_id==20) {
        out << "P\t" << to_string(path_id) << "\t14+,17+\t*\n";
        supported_records.emplace_back("del4");
        supported_records.emplace_back("del5");
    }
    if (path_id==21) {
        out << "P\t" << to_string(path_id) << "\t13+,17+\t*\n";
        supported_records.emplace_back("del3");
        supported_records.emplace_back("del4");
        supported_records.emplace_back("del5");
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
    left_rep_id="rep5";
    right_rep_id="rep7";
    if (path_id==22) {
        out << "P\t" << to_string(path_id) << "\t17+," << left_rep_id << "+,20+\t*\n";
        supported_records.emplace_back(left_rep_id);
        supported_records.emplace_back("del6");
    }
    if (path_id==23) {
        out << "P\t" << to_string(path_id) << "\t18+," << right_rep_id << "+,21+\t*\n";
        supported_records.emplace_back("del6");
        supported_records.emplace_back(right_rep_id);
    }
    if (path_id==24) {
        out << "P\t" << to_string(path_id) << "\t17+," << left_rep_id << "+," << right_rep_id << "+,21+\t*\n";
        supported_records.emplace_back(left_rep_id);
        supported_records.emplace_back("del6");
        supported_records.emplace_back(right_rep_id);
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
    if (path_id==25) {
        out << "P\t" << to_string(path_id) << "\t25+,26-,28+,29+\t*\n";
        supported_records.emplace_back("inv3");
        supported_records.emplace_back("del8");
    }
    if (path_id==26) {
        out << "P\t" << to_string(path_id) << "\t25+,26+,28-,29+\t*\n";
        supported_records.emplace_back("del8");
        supported_records.emplace_back("inv4");
    }
    if (path_id==27) {
        out << "P\t" << to_string(path_id) << "\t25+,26-,28-,29+\t*\n";
        supported_records.emplace_back("inv3");
        supported_records.emplace_back("del8");
        supported_records.emplace_back("inv4");
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
    left_rep_id="rep9";
    if (path_id==28) {
        out << "P\t" << to_string(path_id) << "\t29+,30-," << left_rep_id << "+,32+\t*\n";
        supported_records.emplace_back("inv5");
        supported_records.emplace_back(left_rep_id);
    }
    if (path_id==29) {
        out << "P\t" << to_string(path_id) << "\t30+," << left_rep_id << "+,32-,33+\t*\n";
        supported_records.emplace_back(left_rep_id);
        supported_records.emplace_back("inv6");
    }
    if (path_id==30) {
        out << "P\t" << to_string(path_id) << "\t29+,30-," << left_rep_id << "+,32-,33+\t*\n";
        supported_records.emplace_back("inv5");
        supported_records.emplace_back(left_rep_id);
        supported_records.emplace_back("inv6");
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
    if (path_id==31) {
        out << "P\t" << to_string(path_id) << "\t45+,46-,47-,48+\t*\n";
        supported_records.emplace_back("inv9");
        supported_records.emplace_back("inv10");
    }
    if (path_id==32) {
        out << "P\t" << to_string(path_id) << "\t46+,47-,48-,49+\t*\n";
        supported_records.emplace_back("inv10");
        supported_records.emplace_back("inv11");
    }
    if (path_id==33) {
        out << "P\t" << to_string(path_id) << "\t45+,46-,47-,48-,49+\t*\n";
        supported_records.emplace_back("inv9");
        supported_records.emplace_back("inv10");
        supported_records.emplace_back("inv11");
    }
    if (path_id==34) {
        out << "P\t" << to_string(path_id) << "\t45+,46-,49+\t*\n";
        supported_records.emplace_back("inv9");
        supported_records.emplace_back("inv10");
        supported_records.emplace_back("inv11");
    }
    if (path_id==35) {
        out << "P\t" << to_string(path_id) << "\t45+,48-,49+\t*\n";
        supported_records.emplace_back("inv9");
        supported_records.emplace_back("inv10");
        supported_records.emplace_back("inv11");
    }
}


void print_gfa_colors(ofstream& out) {
    const string COLOR_REF = "gray";
    const string COLOR_INS = "red";
    const string COLOR_REP = "green";

    int32_t i;

    out << "Name,Colour\n";
    for (i=1; i<=49; i++) out << to_string(i) << "," << COLOR_REF << "\n";
    for (i=1; i<=8; i++) out << "ins" << to_string(i) << "," << COLOR_INS << "\n";
    for (i=1; i<=12; i++) out << "rep" << to_string(i) << "," << COLOR_REP << "\n";
}


unordered_map<string,string> get_chromosomes() {
    unordered_map<string,string> out;
    out.emplace("chr1","AAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAA");
    return out;
}


/**
 * Zero-based, non-overlapping, [x..y).
 */
unordered_map<string,vector<interval_t>> get_tandem_track() {
    unordered_map<string,vector<interval_t>> out;
    return out;
}







//void get_edge_record_map(const HashGraph& graph, const vector<string>& node_ids, vector<pair<edge_t,size_t>>& out) {
//    handle_t handle_from, handle_to;
//    out.clear();
//
//    // dup4
//    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1"))+1);
//    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"dup4"))+1);
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),0);  // dup4
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),1);  // dup4_prime
//    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"dup4"))+1);
//    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"2"))+1);
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),0);  // dup4
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),1);  // dup4_prime
//
//    // dup2
//    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"2"))+1);
//    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"dup2"))+1);
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),2);  // dup2
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),3);  // dup2_prime
//    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"dup2"))+1);
//    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"3"))+1);
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),2);  // dup2
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),3);  // dup2_prime
//
//    // dup1
//    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"3"))+1);
//    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"dup1"))+1);
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),4);  // dup1
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),5);  // dup1_prime
//    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"dup1"))+1);
//    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"4"))+1);
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),4);  // dup1
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),5);  // dup1_prime
//
//    // inv5
//    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"4"))+1);
//    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"inv5"))+1);
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),6);  // inv5
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),7);  // inv5_prime
//    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"inv5"))+1);
//    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"11"))+1);
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),6);  // inv5
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),7);  // inv5_prime
//
//    // dup3
//    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"5"))+1);
//    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"dup3"))+1);
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),8);  // dup3
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),9);  // dup3_prime
//    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"dup3"))+1);
//    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"6"))+1);
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),8);  // dup3
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),9);  // dup3_prime
//
//    // inv4
//    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"6"))+1);
//    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"inv4"))+1);
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),10);  // inv4
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),11);  // inv4_prime
//    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"inv4"))+1);
//    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"15"))+1);
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),10);  // inv4
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),11);  // inv4_prime
//
//    // inv2
//    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"7"))+1);
//    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"inv2"))+1);
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),12);  // inv2
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),13);  // inv2_prime
//    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"inv2"))+1);
//    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"10"))+1);
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),12);  // inv2
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),13);  // inv2_prime
//
//    // inv1
//    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"8"))+1);
//    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"inv1"))+1);
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),14);  // inv1
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),15);  // inv1_prime
//    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"inv1"))+1);
//    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"13"))+1);
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),14);  // inv1
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),15);  // inv1_prime
//
//    // inv3
//    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"11"))+1);
//    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"inv3"))+1);
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),16);  // inv3
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),17);  // inv3_prime
//    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"inv3"))+1);
//    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"14"))+1);
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),16);  // inv3
//    out.emplace_back(graph.edge_handle(handle_from,handle_to),17);  // inv3_prime
//}


//void test_vcf_records_with_edges_impl(const string& from, bool from_is_forward, const string& to, bool to_is_forward, vector<edge_t>& edges, vector<VcfRecord>& records, VariantGraph& graph, const vector<string>& node_labels) {
//    handle_t handle_from = graph.graph.get_handle(distance(node_labels.begin(),lower_bound(node_labels.begin(),node_labels.end(),from))+1);
//    if (!from_is_forward) handle_from=graph.graph.flip(handle_from);
//    handle_t handle_to = graph.graph.get_handle(distance(node_labels.begin(),lower_bound(node_labels.begin(),node_labels.end(),to))+1);
//    if (!to_is_forward) handle_to=graph.graph.flip(handle_to);
//    edges.emplace_back(handle_from,handle_to);
//}


//void test_vcf_records_with_edges(VariantGraph& graph, const vector<string>& node_labels) {
//    vector<edge_t> edges;
//    vector<VcfRecord> records;
//    string id;
//
//    id="dup1";
//    edges.clear();
//    test_vcf_records_with_edges_impl("3",true,"dup1",true,edges,records,graph,node_labels);
//    test_vcf_records_with_edges_impl("dup1",true,"4",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=2 || !records.at(0).id.starts_with(id)) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="dup2";
//    edges.clear();
//    test_vcf_records_with_edges_impl("2",true,"dup2",true,edges,records,graph,node_labels);
//    test_vcf_records_with_edges_impl("dup2",true,"3",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=2 || !records.at(0).id.starts_with(id)) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="dup3";
//    edges.clear();
//    test_vcf_records_with_edges_impl("5",true,"dup3",true,edges,records,graph,node_labels);
//    test_vcf_records_with_edges_impl("dup3",true,"6",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=2 || !records.at(0).id.starts_with(id)) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="dup4";
//    edges.clear();
//    test_vcf_records_with_edges_impl("1",true,"dup4",true,edges,records,graph,node_labels);
//    test_vcf_records_with_edges_impl("dup4",true,"2",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=2 || !records.at(0).id.starts_with(id)) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="inv1";
//    edges.clear();
//    test_vcf_records_with_edges_impl("8",true,"inv1",true,edges,records,graph,node_labels);
//    test_vcf_records_with_edges_impl("inv1",true,"13",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=2 || !records.at(0).id.starts_with(id)) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="inv2";
//    edges.clear();
//    test_vcf_records_with_edges_impl("7",true,"inv2",true,edges,records,graph,node_labels);
//    test_vcf_records_with_edges_impl("inv2",true,"10",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=2 || !records.at(0).id.starts_with(id)) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="inv3";
//    edges.clear();
//    test_vcf_records_with_edges_impl("11",true,"inv3",true,edges,records,graph,node_labels);
//    test_vcf_records_with_edges_impl("inv3",true,"14",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=2 || !records.at(0).id.starts_with(id)) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="inv4";
//    edges.clear();
//    test_vcf_records_with_edges_impl("6",true,"inv4",true,edges,records,graph,node_labels);
//    test_vcf_records_with_edges_impl("inv4",true,"15",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=2 || !records.at(0).id.starts_with(id)) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="inv5";
//    edges.clear();
//    test_vcf_records_with_edges_impl("4",true,"inv5",true,edges,records,graph,node_labels);
//    test_vcf_records_with_edges_impl("inv5",true,"11",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=2 || !records.at(0).id.starts_with(id)) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//}


int main(int argc, char* argv[]) {
    const path ROOT_DIR = path(argv[1]);
    const int32_t CLOSURE_1 = stoi(argv[2]);
    const int32_t CLOSURE_2 = stoi(argv[3]);

    const path INPUT_VCF = ROOT_DIR/"input.vcf";
    const path TRUTH_GFA = ROOT_DIR/"truth.gfa";
    const path TRUTH_GFA_COLORS = ROOT_DIR/"truth.csv";
    const path TEST_GFA = ROOT_DIR/"test.gfa";
    const int32_t SIGNATURE_N_STEPS = 10;
    const int32_t FLANK_LENGTH = 10;
    const int32_t INTERIOR_FLANK_LENGTH = 10;

    const unordered_map<string,string> chromosomes = get_chromosomes();
    const unordered_map<string,vector<interval_t>> tandem_track = get_tandem_track();
    size_t i, j;
    size_t n_records;
    string command;
    vector<string> supported_records;
    vector<VcfRecord> records;

    // Printing truth files
    ofstream input_vcf(INPUT_VCF.string());
    print_truth_vcf_header(input_vcf);
    print_truth_vcf(input_vcf);
    input_vcf.close();
    ofstream truth_gfa(TRUTH_GFA.string());
    print_truth_gfa(truth_gfa,CLOSURE_1==1,CLOSURE_2==1,0,supported_records);
    truth_gfa.close();
    ofstream truth_gfa_colors(TRUTH_GFA_COLORS.string());
    print_gfa_colors(truth_gfa_colors);
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
    for (i=1; i<=35; i++) {
        ofstream tmp_gfa("tmp.gfa");
        print_truth_gfa(tmp_gfa,true,true,i,supported_records);
        tmp_gfa.close();
        graph.load_gfa("tmp.gfa");
        ofstream supported_test_vcf("supported_test.vcf"); ofstream unsupported_test_vcf("unsupported_test.vcf");
        graph.print_supported_vcf_records(supported_test_vcf,unsupported_test_vcf,false);
        supported_test_vcf.close(); unsupported_test_vcf.close();
        ofstream supported_truth_vcf("supported_truth.vcf");
        for (auto& id: supported_records) {
            for (j=0; j<n_records; j++) {
                if (records.at(j).id==id) {
                    records.at(j).print(supported_truth_vcf);
                    supported_truth_vcf << '\n';
                    break;
                }
            }
        }
        supported_truth_vcf.close();
        command.clear(); command.append("sort supported_test.vcf > supported_test_sorted.vcf"); run_command(command);
        command.clear(); command.append("sort supported_truth.vcf > supported_truth_sorted.vcf"); run_command(command);
        command.clear(); command.append("diff --brief supported_test_sorted.vcf supported_truth_sorted.vcf"); run_command(command);
    }

    cerr << "Removing temporary files...\n";
    command.clear(); command.append("rm -f tmp*.txt tmp*.gfa supported_*.vcf unsupported_*.vcf"); run_command(command);
}