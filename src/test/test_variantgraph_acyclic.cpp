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

    out << "chr1\t19\tdup4\tA\t<DUP>\t" << INFIX << "\tSVTYPE=DUP;SVLEN=71;\t" << SUFFIX;
    out << "chr1\t19\tdup4_prime\tA\t<CNV>\t" << INFIX << "\tSVTYPE=CNV;SVLEN=71;\t" << SUFFIX;
    out << "chr1\t29\tdup2\tA\t<DUP>\t" << INFIX << "\tSVTYPE=DUP;SVLEN=21;\t" << SUFFIX;
    out << "chr1\t29\tdup2_prime\tA\t<CNV>\t" << INFIX << "\tSVTYPE=CNV;SVLEN=21;\t" << SUFFIX;
    out << "chr1\t39\tdup1\tA\t<DUP>\t" << INFIX << "\tSVTYPE=DUP;SVLEN=31;\t" << SUFFIX;
    out << "chr1\t39\tdup1_prime\tA\t<CNV>\t" << INFIX << "\tSVTYPE=CNV;SVLEN=31;\t" << SUFFIX;
    out << "chr1\t51\tinv5\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=84;\t" << SUFFIX;
    out << "chr1\t51\tinv5_prime\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=84;\t" << SUFFIX;
    out << "chr1\t54\tdup3\tA\t<DUP>\t" << INFIX << "\tSVTYPE=DUP;SVLEN=26;\t" << SUFFIX;
    out << "chr1\t54\tdup3_prime\tA\t<CNV>\t" << INFIX << "\tSVTYPE=CNV;SVLEN=26;\t" << SUFFIX;
    out << "chr1\t99\tinv4\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=71;\t" << SUFFIX;
    out << "chr1\t99\tinv4_prime\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=71;\t" << SUFFIX;
    out << "chr1\t109\tinv2\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=21;\t" << SUFFIX;
    out << "chr1\t109\tinv2_prime\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=21;\t" << SUFFIX;
    out << "chr1\t119\tinv1\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=31;\t" << SUFFIX;
    out << "chr1\t119\tinv1_prime\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=31;\t" << SUFFIX;
    out << "chr1\t139\tinv3\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=21;\t" << SUFFIX;
    out << "chr1\t139\tinv3_prime\tA\t<INV>\t" << INFIX << "\tSVTYPE=INV;SVLEN=21;\t" << SUFFIX;
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


void get_edge_record_map(const HashGraph& graph, const vector<string>& node_ids, vector<pair<edge_t,size_t>>& out) {
    handle_t handle_from, handle_to;
    out.clear();

    // dup4
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"1"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"dup4"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),0);  // dup4
    out.emplace_back(graph.edge_handle(handle_from,handle_to),1);  // dup4_prime
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"dup4"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"2"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),0);  // dup4
    out.emplace_back(graph.edge_handle(handle_from,handle_to),1);  // dup4_prime

    // dup2
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"2"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"dup2"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),2);  // dup2
    out.emplace_back(graph.edge_handle(handle_from,handle_to),3);  // dup2_prime
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"dup2"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"3"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),2);  // dup2
    out.emplace_back(graph.edge_handle(handle_from,handle_to),3);  // dup2_prime

    // dup1
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"3"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"dup1"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),4);  // dup1
    out.emplace_back(graph.edge_handle(handle_from,handle_to),5);  // dup1_prime
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"dup1"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"4"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),4);  // dup1
    out.emplace_back(graph.edge_handle(handle_from,handle_to),5);  // dup1_prime

    // inv5
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"4"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"inv5"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),6);  // inv5
    out.emplace_back(graph.edge_handle(handle_from,handle_to),7);  // inv5_prime
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"inv5"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"11"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),6);  // inv5
    out.emplace_back(graph.edge_handle(handle_from,handle_to),7);  // inv5_prime

    // dup3
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"5"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"dup3"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),8);  // dup3
    out.emplace_back(graph.edge_handle(handle_from,handle_to),9);  // dup3_prime
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"dup3"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"6"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),8);  // dup3
    out.emplace_back(graph.edge_handle(handle_from,handle_to),9);  // dup3_prime

    // inv4
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"6"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"inv4"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),10);  // inv4
    out.emplace_back(graph.edge_handle(handle_from,handle_to),11);  // inv4_prime
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"inv4"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"15"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),10);  // inv4
    out.emplace_back(graph.edge_handle(handle_from,handle_to),11);  // inv4_prime

    // inv2
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"7"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"inv2"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),12);  // inv2
    out.emplace_back(graph.edge_handle(handle_from,handle_to),13);  // inv2_prime
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"inv2"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"10"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),12);  // inv2
    out.emplace_back(graph.edge_handle(handle_from,handle_to),13);  // inv2_prime

    // inv1
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"8"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"inv1"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),14);  // inv1
    out.emplace_back(graph.edge_handle(handle_from,handle_to),15);  // inv1_prime
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"inv1"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"13"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),14);  // inv1
    out.emplace_back(graph.edge_handle(handle_from,handle_to),15);  // inv1_prime

    // inv3
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"11"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"inv3"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),16);  // inv3
    out.emplace_back(graph.edge_handle(handle_from,handle_to),17);  // inv3_prime
    handle_from=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"inv3"))+1);
    handle_to=graph.get_handle(distance(node_ids.begin(),lower_bound(node_ids.begin(),node_ids.end(),"14"))+1);
    out.emplace_back(graph.edge_handle(handle_from,handle_to),16);  // inv3
    out.emplace_back(graph.edge_handle(handle_from,handle_to),17);  // inv3_prime
}


void test_vcf_records_with_edges_impl(const string& from, bool from_is_forward, const string& to, bool to_is_forward, vector<edge_t>& edges, vector<VcfRecord>& records, VariantGraph& graph, const vector<string>& node_labels) {
    handle_t handle_from = graph.graph.get_handle(distance(node_labels.begin(),lower_bound(node_labels.begin(),node_labels.end(),from))+1);
    if (!from_is_forward) handle_from=graph.graph.flip(handle_from);
    handle_t handle_to = graph.graph.get_handle(distance(node_labels.begin(),lower_bound(node_labels.begin(),node_labels.end(),to))+1);
    if (!to_is_forward) handle_to=graph.graph.flip(handle_to);
    edges.emplace_back(handle_from,handle_to);
}


void test_vcf_records_with_edges(VariantGraph& graph, const vector<string>& node_labels) {
    vector<edge_t> edges;
    vector<VcfRecord> records;
    string id;

    id="dup1";
    edges.clear();
    test_vcf_records_with_edges_impl("3",true,"dup1",true,edges,records,graph,node_labels);
    test_vcf_records_with_edges_impl("dup1",true,"4",true,edges,records,graph,node_labels);
    graph.get_vcf_records_with_edges(edges,records);
    if (records.size()!=2 || !records.at(0).id.starts_with(id)) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);

    id="dup2";
    edges.clear();
    test_vcf_records_with_edges_impl("2",true,"dup2",true,edges,records,graph,node_labels);
    test_vcf_records_with_edges_impl("dup2",true,"3",true,edges,records,graph,node_labels);
    graph.get_vcf_records_with_edges(edges,records);
    if (records.size()!=2 || !records.at(0).id.starts_with(id)) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);

    id="dup3";
    edges.clear();
    test_vcf_records_with_edges_impl("5",true,"dup3",true,edges,records,graph,node_labels);
    test_vcf_records_with_edges_impl("dup3",true,"6",true,edges,records,graph,node_labels);
    graph.get_vcf_records_with_edges(edges,records);
    if (records.size()!=2 || !records.at(0).id.starts_with(id)) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);

    id="dup4";
    edges.clear();
    test_vcf_records_with_edges_impl("1",true,"dup4",true,edges,records,graph,node_labels);
    test_vcf_records_with_edges_impl("dup4",true,"2",true,edges,records,graph,node_labels);
    graph.get_vcf_records_with_edges(edges,records);
    if (records.size()!=2 || !records.at(0).id.starts_with(id)) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);

    id="inv1";
    edges.clear();
    test_vcf_records_with_edges_impl("8",true,"inv1",true,edges,records,graph,node_labels);
    test_vcf_records_with_edges_impl("inv1",true,"13",true,edges,records,graph,node_labels);
    graph.get_vcf_records_with_edges(edges,records);
    if (records.size()!=2 || !records.at(0).id.starts_with(id)) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);

    id="inv2";
    edges.clear();
    test_vcf_records_with_edges_impl("7",true,"inv2",true,edges,records,graph,node_labels);
    test_vcf_records_with_edges_impl("inv2",true,"10",true,edges,records,graph,node_labels);
    graph.get_vcf_records_with_edges(edges,records);
    if (records.size()!=2 || !records.at(0).id.starts_with(id)) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);

    id="inv3";
    edges.clear();
    test_vcf_records_with_edges_impl("11",true,"inv3",true,edges,records,graph,node_labels);
    test_vcf_records_with_edges_impl("inv3",true,"14",true,edges,records,graph,node_labels);
    graph.get_vcf_records_with_edges(edges,records);
    if (records.size()!=2 || !records.at(0).id.starts_with(id)) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);

    id="inv4";
    edges.clear();
    test_vcf_records_with_edges_impl("6",true,"inv4",true,edges,records,graph,node_labels);
    test_vcf_records_with_edges_impl("inv4",true,"15",true,edges,records,graph,node_labels);
    graph.get_vcf_records_with_edges(edges,records);
    if (records.size()!=2 || !records.at(0).id.starts_with(id)) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);

    id="inv5";
    edges.clear();
    test_vcf_records_with_edges_impl("4",true,"inv5",true,edges,records,graph,node_labels);
    test_vcf_records_with_edges_impl("inv5",true,"11",true,edges,records,graph,node_labels);
    graph.get_vcf_records_with_edges(edges,records);
    if (records.size()!=2 || !records.at(0).id.starts_with(id)) throw runtime_error("get_vcf_records_with_edges() failed on VCF record "+id);
}


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
    const size_t n_records = records.size();
    VariantGraph graph(chromosomes,tandem_track);

    cerr << "Testing acyclic GFA...\n";
    graph.build(records,FLANK_LENGTH,INTERIOR_FLANK_LENGTH,INT32_MAX,INT32_MAX,false,{},true);
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

    cerr << "Testing edge-record map (1/2)...\n";
    unordered_map<edge_t,vector<size_t>> map1 = graph.get_edge_record_map();
    vector<string> node_labels=graph.load_gfa(TRUTH_GFA.string());
    vector<pair<edge_t,size_t>> edge_record_map;
    get_edge_record_map(graph.graph,node_labels,edge_record_map);
    graph.load_edge_record_map(edge_record_map,n_records);
    test_vcf_records_with_edges(graph,node_labels);

    cerr << "Testing edge-record map (2/2)...\n";
    ofstream map1_txt("tmp1.txt");
    for (const auto& pair: map1) {
        for (const auto& id: pair.second) { map1_txt << to_string(id) << ','; }
        map1_txt << '\n';
    }
    map1_txt.close();
    unordered_map<edge_t,vector<size_t>> map2 = graph.get_edge_record_map();
    ofstream map2_txt("tmp2.txt");
    for (const auto& pair: map2) {
        for (const auto& id: pair.second) { map2_txt << to_string(id) << ','; }
        map2_txt << '\n';
    }
    map2_txt.close();
    command.clear(); command.append("sort tmp1.txt > tmp1_sorted.txt"); run_command(command);
    command.clear(); command.append("sort tmp2.txt > tmp2_sorted.txt"); run_command(command);
    command.clear(); command.append("diff --brief tmp1_sorted.txt tmp2_sorted.txt"); run_command(command);

    cerr << "Removing temporary files...\n";
    command.clear(); command.append("rm -f input*.vcf truth*.gfa test*.gfa tmp*.txt"); run_command(command);
}