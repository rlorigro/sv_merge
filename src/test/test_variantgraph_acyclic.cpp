#include "VcfReader.hpp"
#include "VariantGraph.hpp"
#include "misc.hpp"

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

    out << "chr1\t20\tdup4\tA\t<DUP>\t" << INFIX << "\tSVTYPE=DUP;SVLEN=71;\t" << SUFFIX;
    out << "chr1\t20\tdup4_prime\tA\t<CNV>\t" << INFIX << "\tSVTYPE=CNV;SVLEN=71;\t" << SUFFIX;
    out << "chr1\t30\tdup2\tA\t<DUP>\t" << INFIX << "\tSVTYPE=DUP;SVLEN=21;\t" << SUFFIX;
    out << "chr1\t30\tdup2_prime\tA\t<CNV>\t" << INFIX << "\tSVTYPE=CNV;SVLEN=21;\t" << SUFFIX;
    out << "chr1\t40\tdup1\tA\t<DUP>\t" << INFIX << "\tSVTYPE=DUP;SVLEN=31;\t" << SUFFIX;
    out << "chr1\t40\tdup1_prime\tA\t<CNV>\t" << INFIX << "\tSVTYPE=CNV;SVLEN=31;\t" << SUFFIX;
    out << "chr1\t52\tinv5\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=84;\t" << SUFFIX;
    out << "chr1\t52\tinv5_prime\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=84;\t" << SUFFIX;
    out << "chr1\t55\tdup3\tA\t<DUP>\t" << INFIX << "\tSVTYPE=DUP;SVLEN=26;\t" << SUFFIX;
    out << "chr1\t55\tdup3_prime\tA\t<CNV>\t" << INFIX << "\tSVTYPE=CNV;SVLEN=26;\t" << SUFFIX;
    out << "chr1\t100\tinv4\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=71;\t" << SUFFIX;
    out << "chr1\t100\tinv4_prime\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=71;\t" << SUFFIX;
    out << "chr1\t110\tinv2\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=21;\t" << SUFFIX;
    out << "chr1\t110\tinv2_prime\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=21;\t" << SUFFIX;
    out << "chr1\t120\tinv1\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=31;\t" << SUFFIX;
    out << "chr1\t120\tinv1_prime\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=31;\t" << SUFFIX;
    out << "chr1\t140\tinv3\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=21;\t" << SUFFIX;
    out << "chr1\t140\tinv3_prime\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=21;\t" << SUFFIX;
}


void print_truth_gfa(ofstream& out) {
    // Ref
    out << "S\t1\tAAAAAAAAAAAAAAAAAAA\n";
    out << "S\t2\tAAAAAAAAAA\n";
    out << "S\t3\tAAAAAAAAAA\n";
    out << "S\t4\tAAAAAAAAAAAA\n";
    out << "S\t5\tAAA\n";
    out << "S\t6\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";
    out << "S\t7\tAAAAAAAAAA\n";
    out << "S\t8\tAAAAAAAAAA\n";
    out << "S\t9\tAAAAAAAAAAA\n";
    out << "S\t10\tAAAAA\n";
    out << "S\t11\tAAAA\n";
    out << "S\t12\tAAAAAAAAAAA\n";
    out << "S\t13\tAAAAAAAAAA\n";
    out << "S\t14\tAAAAAAAAAA\n";
    out << "S\t15\tAAAAA\n";

    // Non-ref
    out << "S\tdup1\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";
    out << "S\tdup2\tAAAAAAAAAAAAAAAAAAAAA\n";
    out << "S\tdup3\tAAAAAAAAAAAAAAAAAAAAAAAAAA\n";
    out << "S\tdup4\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";
    out << "S\tinv5\tTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n";
    out << "S\tinv4\tTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n";
    out << "S\tinv2\tTTTTTTTTTTTTTTTTTTTTT\n";
    out << "S\tinv3\tTTTTTTTTTTTTTTTTTTTTT\n";
    out << "S\tinv1\tTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n";

    // Ref
    out << "L\t1\t+\t2\t+\t*\n";
    out << "L\t2\t+\t3\t+\t*\n";
    out << "L\t3\t+\t4\t+\t*\n";
    out << "L\t4\t+\t5\t+\t*\n";
    out << "L\t5\t+\t6\t+\t*\n";
    out << "L\t6\t+\t7\t+\t*\n";
    out << "L\t7\t+\t8\t+\t*\n";
    out << "L\t8\t+\t9\t+\t*\n";
    out << "L\t9\t+\t10\t+\t*\n";
    out << "L\t10\t+\t11\t+\t*\n";
    out << "L\t11\t+\t12\t+\t*\n";
    out << "L\t12\t+\t13\t+\t*\n";
    out << "L\t13\t+\t14\t+\t*\n";
    out << "L\t14\t+\t15\t+\t*\n";

    // Non-ref
    out << "L\t1\t+\tdup4\t+\t*\n";
    out << "L\tdup4\t+\t2\t+\t*\n";
    out << "L\t2\t+\tdup2\t+\t*\n";
    out << "L\tdup2\t+\t3\t+\t*\n";
    out << "L\t3\t+\tdup1\t+\t*\n";
    out << "L\tdup1\t+\t4\t+\t*\n";
    out << "L\t4\t+\tinv5\t+\t*\n";
    out << "L\t5\t+\tdup3\t+\t*\n";
    out << "L\tdup3\t+\t6\t+\t*\n";
    out << "L\t6\t+\tinv4\t+\t*\n";
    out << "L\t7\t+\tinv2\t+\t*\n";
    out << "L\t8\t+\tinv1\t+\t*\n";
    out << "L\tinv2\t+\t10\t+\t*\n";
    out << "L\tinv5\t+\t11\t+\t*\n";
    out << "L\t11\t+\tinv3\t+\t*\n";
    out << "L\tinv1\t+\t13\t+\t*\n";
    out << "L\tinv3\t+\t14\t+\t*\n";
    out << "L\tinv4\t+\t15\t+\t*\n";
}


unordered_map<string,string> get_chromosomes() {
    unordered_map<string,string> out;
    out.emplace("chr1","AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    return out;
}


/**
 * Zero-based, non-overlapping, [x..y).
 */
unordered_map<string,vector<interval_t>> get_tandem_track() {
    unordered_map<string,vector<interval_t>> out;
    return out;
}




//
//void test_vcf_records_with_edges_impl(const string& from, bool from_is_forward, const string& to, bool to_is_forward, vector<edge_t>& edges, vector<VcfRecord>& records, VariantGraph& graph, const vector<string>& node_labels) {
//    handle_t handle_from = graph.graph.get_handle(distance(node_labels.begin(),lower_bound(node_labels.begin(),node_labels.end(),from))+1);
//    if (!from_is_forward) handle_from=graph.graph.flip(handle_from);
//    handle_t handle_to = graph.graph.get_handle(distance(node_labels.begin(),lower_bound(node_labels.begin(),node_labels.end(),to))+1);
//    if (!to_is_forward) handle_to=graph.graph.flip(handle_to);
//    edges.emplace_back(handle_from,handle_to);
//}
//
//
//void test_vcf_records_with_edges(VariantGraph& graph, const vector<string>& node_labels) {
//    vector<edge_t> edges;
//    vector<VcfRecord> records;
//    string id;
//
//    id="sniffles_del1";
//    edges.clear();
//    test_vcf_records_with_edges_impl("1_2",true,"1_9",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="pbsv_del2";
//    edges.clear();
//    test_vcf_records_with_edges_impl("1_4",true,"1_13",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="pbsv_del3";
//    edges.clear();
//    test_vcf_records_with_edges_impl("1_7",true,"1_17",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="sniffles_inv1";
//    edges.clear();
//    test_vcf_records_with_edges_impl("1_3",true,"1_18",false,edges,records,graph,node_labels);
//    test_vcf_records_with_edges_impl("1_4",false,"1_19",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="pbsv_inv2";
//    edges.clear();
//    test_vcf_records_with_edges_impl("1_6",true,"1_21",false,edges,records,graph,node_labels);
//    test_vcf_records_with_edges_impl("1_7",false,"1_22",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="pav_inv3";
//    edges.clear();
//    test_vcf_records_with_edges_impl("1_10",true,"1_24",false,edges,records,graph,node_labels);
//    test_vcf_records_with_edges_impl("1_11",false,"1_25",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="svim_inv4";
//    edges.clear();
//    test_vcf_records_with_edges_impl("1_20",true,"1_28",false,edges,records,graph,node_labels);
//    test_vcf_records_with_edges_impl("1_21",false,"1_29",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="pbsv_ins1";
//    edges.clear();
//    test_vcf_records_with_edges_impl("1_4",true,"ins1",true,edges,records,graph,node_labels);
//    test_vcf_records_with_edges_impl("ins1",true,"1_6",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="pav_ins2";
//    edges.clear();
//    test_vcf_records_with_edges_impl("1_9",true,"ins2",true,edges,records,graph,node_labels);
//    test_vcf_records_with_edges_impl("ins2",true,"1_10",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="pav_ins3";
//    edges.clear();
//    test_vcf_records_with_edges_impl("1_11",true,"ins3",true,edges,records,graph,node_labels);
//    test_vcf_records_with_edges_impl("ins3",true,"1_12",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="pbsv_rep1";
//    edges.clear();
//    test_vcf_records_with_edges_impl("1_8",true,"rep1",true,edges,records,graph,node_labels);
//    test_vcf_records_with_edges_impl("rep1",true,"1_23",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    edges.clear();
//    test_vcf_records_with_edges_impl("1_26",true,"1_2",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=3) throw runtime_error("get_vcf_records_with_edges() failed (1) on VCF records sniffles_dup1 = sniffles_bnd19");
//    vector<string> ids;
//    for (const auto& record: records) ids.emplace_back(record.id);
//    std::sort(ids.begin(),ids.end());
//    if (ids.at(0)!="cutesv_bnd20" || ids.at(1)!="sniffles_bnd19" || ids.at(2)!="sniffles_dup1") throw runtime_error("get_vcf_records_with_edges() failed (2) on VCF records sniffles_dup1 = sniffles_bnd19");
//
//    id="svim_bnd10";
//    edges.clear();
//    test_vcf_records_with_edges_impl("bnd10",true,"1_18",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="svim_bnd11";
//    edges.clear();
//    test_vcf_records_with_edges_impl("1_25",true,"bnd11",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="cutesv_bnd12";
//    edges.clear();
//    test_vcf_records_with_edges_impl("bnd12",true,"1_26",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="sniffles_bnd1";
//    edges.clear();
//    test_vcf_records_with_edges_impl("2_40",true,"1_3",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="pbsv_bnd2";
//    edges.clear();
//    test_vcf_records_with_edges_impl("1_7",true,"2_50",false,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="pav_bnd3";
//    edges.clear();
//    test_vcf_records_with_edges_impl("1_8",false,"2_58",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="svim_bnd4";
//    edges.clear();
//    test_vcf_records_with_edges_impl("1_19",true,"2_61",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="svim_bnd5";
//    edges.clear();
//    test_vcf_records_with_edges_impl("bnd5",true,"1_20",true,edges,records,graph,node_labels);
//    test_vcf_records_with_edges_impl("2_172",false,"bnd5",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="svim_bnd6";
//    edges.clear();
//    test_vcf_records_with_edges_impl("3_67",true,"bnd6",true,edges,records,graph,node_labels);
//    test_vcf_records_with_edges_impl("bnd6",true,"1_24",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="cutesv_bnd7";
//    edges.clear();
//    test_vcf_records_with_edges_impl("1_129",true,"bnd7",true,edges,records,graph,node_labels);
//    test_vcf_records_with_edges_impl("bnd7",true,"3_67",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="cutesv_bnd8";
//    edges.clear();
//    test_vcf_records_with_edges_impl("1_130",true,"bnd8",true,edges,records,graph,node_labels);
//    test_vcf_records_with_edges_impl("bnd8",true,"2_285",false,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="pav_bnd9";
//    edges.clear();
//    test_vcf_records_with_edges_impl("1_13",true,"2_167",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="svision_ins0";
//    edges.clear();
//    test_vcf_records_with_edges_impl("ins0",true,"1_1",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="svision_ins5";
//    edges.clear();
//    test_vcf_records_with_edges_impl("1_130",true,"ins5",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="svision_dup2";
//    edges.clear();
//    test_vcf_records_with_edges_impl("1_12",true,"1_1",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="svision_dup3";
//    edges.clear();
//    test_vcf_records_with_edges_impl("1_130",true,"1_130",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="svision_inv5";
//    edges.clear();
//    test_vcf_records_with_edges_impl("1_1",false,"1_17",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="svision_inv6";
//    edges.clear();
//    test_vcf_records_with_edges_impl("1_129",true,"1_130",false,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="svision_rep2";
//    edges.clear();
//    test_vcf_records_with_edges_impl("rep2",true,"1_18",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="svision_rep3";
//    edges.clear();
//    test_vcf_records_with_edges_impl("1_129",true,"rep3",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//
//    id="svision_bnd21";
//    edges.clear();
//    test_vcf_records_with_edges_impl("2_285",true,"1_1",true,edges,records,graph,node_labels);
//    graph.get_vcf_records_with_edges(edges,records);
//    if (records.size()!=1 || records.at(0).id!=id) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
//}














int main(int argc, char* argv[]) {
    const path ROOT_DIR = path(argv[1]);

    const path INPUT_VCF = ROOT_DIR/"input.vcf";
    const path TRUTH_GFA = ROOT_DIR/"truth.gfa";
    const path TEST_GFA = ROOT_DIR/"test.gfa";
    const int32_t FLANK_LENGTH = INT32_MAX;
    const int32_t INTERIOR_FLANK_LENGTH = INT32_MAX;
    const int32_t SIGNATURE_N_STEPS = 10;

    ofstream input_vcf(INPUT_VCF.string());
    print_truth_vcf_header(input_vcf);
    print_truth_vcf(input_vcf);
    input_vcf.close();
    ofstream truth_gfa(TRUTH_GFA.string());
    print_truth_gfa(truth_gfa);
    truth_gfa.close();

    const unordered_map<string,string> chromosomes = get_chromosomes();
    const unordered_map<string,vector<interval_t>> tandem_track = get_tandem_track();
    string command;

    vector<VcfRecord> records;
    VcfReader reader(INPUT_VCF);
    reader.for_record_in_vcf([&](VcfRecord& record) {
        if ( (record.sv_type==VcfReader::TYPE_INSERTION && record.is_symbolic) ||
             ((record.sv_type==VcfReader::TYPE_DELETION || record.sv_type==VcfReader::TYPE_INVERSION || record.sv_type==VcfReader::TYPE_DUPLICATION || record.sv_type==VcfReader::TYPE_REPLACEMENT) && record.sv_length==INT32_MAX)
           ) return;
        records.push_back(record);
    });
    VariantGraph graph(chromosomes,tandem_track);

    cerr << "Testing acyclic GFA...\n";
    graph.build(records,FLANK_LENGTH,INTERIOR_FLANK_LENGTH,INT32_MAX,INT32_MAX,false);
    ofstream test_gfa(TEST_GFA.string());
    graph.to_gfa(TEST_GFA);
    test_gfa.close();

    cerr << "Testing nodes...\n";
    command.clear(); command.append("grep ^S "+TRUTH_GFA.string()+" | cut -f 1,3 | sort > tmp1.txt"); run_command(command);
    command.clear(); command.append("grep ^S "+TEST_GFA.string()+" | cut -f 1,3 | sort > tmp2.txt"); run_command(command);
    command.clear(); command.append("diff --brief tmp1.txt tmp2.txt"); run_command(command);

    cerr << "Testing n. edges...\n";
    command.clear(); command.append("grep ^L "+TRUTH_GFA.string()+" | wc -l > tmp1.txt"); run_command(command);
    command.clear(); command.append("grep ^L "+TEST_GFA.string()+" | wc -l > tmp2.txt"); run_command(command);
    command.clear(); command.append("diff --brief tmp1.txt tmp2.txt"); run_command(command);

    cerr << "Testing local topology (" << SIGNATURE_N_STEPS << " steps)...\n";
    graph.print_graph_signature(SIGNATURE_N_STEPS,"tmp1.txt");
    graph.load_gfa(TRUTH_GFA);
    graph.print_graph_signature(SIGNATURE_N_STEPS,"tmp2.txt");
    command.clear(); command.append("diff --brief tmp1.txt tmp2.txt"); run_command(command);













//    cerr << "Testing vcf_records_with_edges()...\n";
//    truth_vcf.clear(); truth_vcf.open(TRUTH_VCF.string());
//    print_truth_vcf_header(truth_vcf);
//    print_truth_vcf(0,false/*Arbitrary*/,truth_vcf);
//    truth_vcf.close();
//    records.clear();
//    VcfReader reader3(TRUTH_VCF);
//    reader3.for_record_in_vcf([&](VcfRecord& record) {
//        if ( (record.sv_type==VcfReader::TYPE_INSERTION && record.is_symbolic) ||
//             ((record.sv_type==VcfReader::TYPE_DELETION || record.sv_type==VcfReader::TYPE_INVERSION || record.sv_type==VcfReader::TYPE_DUPLICATION || record.sv_type==VcfReader::TYPE_REPLACEMENT) && record.sv_length==INT32_MAX)
//           ) return;
//        records.push_back(record);
//    });
//    graph.build(records,FLANK_LENGTH,INTERIOR_FLANK_LENGTH,INT32_MAX,INT32_MAX,false,caller_ids);
//    truth_gfa.clear(); truth_gfa.open(TRUTH_GFA.string());
//    print_truth_gfa(true,true,true,true,true,true,truth_gfa);
//    truth_gfa.close();
//    node_labels=graph.load_gfa(TRUTH_GFA.string());
//    edge_record_map=get_edge_record_map(graph.graph,node_labels);
//    graph.load_edge_record_map(edge_record_map,get_n_vcf_records());
//    test_vcf_records_with_edges(graph,node_labels);

    cerr << "Removing temporary files...\n";
    command.clear(); command.append("rm -f input*.vcf truth*.gfa test*.gfa tmp*.txt"); run_command(command);
}