#include "fasta.hpp"
#include "CLI11.hpp"
#include "VcfReader.hpp"

#include <unordered_map>
#include <exception>
#include <iostream>

using std::unordered_map;
using std::exception;
using std::ifstream;
using std::ofstream;
using std::cerr;

using namespace sv_merge;


/**
 * Sets `header` to an updated version of the header of `input_vcf`.
 */
void get_vcf_header(const path& input_vcf, string& header) {
    bool is_last_line;
    char c;
    size_t i;
    string line;
    ifstream file;

    file.open(input_vcf,std::ifstream::in);
    if (!file.is_open() || !file.good()) throw runtime_error("ERROR: could not read file: "+input_vcf.string());
    header.clear();
    while (!file.eof()) {
        line.clear();
        c=(char)file.peek();
        if (c!=VcfReader::VCF_COMMENT) break;
        while (file.get(c)) {
            if (c==VcfReader::LINE_END || c==EOF) break;
            else line.push_back(c);
        }
        is_last_line=line.length()>=VcfReader::VCF_HEADER_PREFIX_LENGTH;
        if (is_last_line) {
            for (i=0; i<VcfReader::VCF_HEADER_PREFIX_LENGTH; i++) {
                if (toupper(line.at(i))!=VcfReader::VCF_HEADER_PREFIX.at(i)) { is_last_line=false; break; }
            }
        }
        if (is_last_line) header.append("##INFO=<ID=MATEID,Number=1,Type=String,Description=\"ID of mate breakends\">\n");
        header.append(line+"\n");
    }
    file.close();
}


/**
 * - Removes symbolic breakends and virtual telomeric breakends.
 * - Ensures that REF and ALT in BND records have the same character at POS, and that such a character is the same as
 *   in the ref (this may not be true in e.g. sniffles).
 * - Ensures that every non-single BND record has a symmetrical mate record.
 * - Discards BND records whose ALT allele does not conform to the VCF spec (e.g. sniffles can output ACNNNNNNNNNNNNNNN
 *   in the ALT of a BND record).
 * - Transforms long non-BND calls into sets of BNDs.
 *
 * Remark: the output VCF is not sorted: if it is used for inter-sample SV merging, it should be sorted. If it is used
 * for inter-sample SV merging or for intra-sample SV filtering, it should be post-processed to make a decision on a
 * large call from information on its BNDs.
 */
int main (int argc, char* argv[]) {
    bool has_mate;
    char c;
    int32_t pos2, n_records, min_sv_length, single_breakend_length, chromosome_length;
    size_t alt_length;
    string header, current_chromosome, chromosome2, buffer;
    path ref_fasta, input_vcf, output_vcf;

    const string MATEID_STR = "MATEID";
    const string NEW_MATE_PREFIX = "mate_of_";
    const string SECOND_ADJACENCY_SUFFIX = "_prime";

    CLI::App app{"App description"};
    app.add_option("--input_vcf",input_vcf,"Path to input VCF file, which might contain both BND and non-BND records.")->required();
    app.add_option("--ref_fasta",ref_fasta,"Path to reference sequence FASTA file that corresponds to VCF")->required();
    app.add_option("--min_sv_length",min_sv_length,"Only calls this length or longer are transformed into BNDs. Every other call is kept intact and printed to the output.")->required();
    app.add_option("--single_breakend_length",single_breakend_length,"A long INS is converted to a pair of single breakends with this number of bases each.")->required();
    app.add_option("--output_vcf",output_vcf,"Path to output VCF file")->required();
    CLI11_PARSE(app,argc,argv)
    if (min_sv_length<(single_breakend_length<<1)) throw runtime_error("The following command-line arguments are inconsistent: min_sv_length="+std::to_string(min_sv_length)+" single_breakend_length="+std::to_string(single_breakend_length));

    cerr << "Loading the reference...\n";
    unordered_map<string,string> ref_sequences;
    for_sequence_in_fasta_file(ref_fasta, [&](const Sequence& s){ ref_sequences[s.name] = s.sequence; });
    cerr << "done\n";

    cerr << "Reading the VCF... " << '\n';
    VcfReader vcf_reader(input_vcf);
    vcf_reader.min_qual = numeric_limits<float>::min();
    vcf_reader.min_sv_length = 0;
    vcf_reader.progress_n_lines = 100'000;
    ofstream outstream(output_vcf.string());
    if (!outstream.good() || !outstream.is_open()) throw runtime_error("ERROR: file could not be written: "+output_vcf.string());
    get_vcf_header(input_vcf,header);
    outstream << header;
    n_records=0;
    vcf_reader.for_record_in_vcf([&](VcfRecord& record) {
        n_records++;
        if (record.sv_type==-1) { record.print(outstream); outstream << "\n"; return; }
        else if (record.sv_type==VcfReader::TYPE_BREAKEND) {
            if (record.is_symbolic || !record.is_valid_bnd_alt() || record.is_breakend_virtual(ref_sequences)) return;
            has_mate=record.get_info_field(MATEID_STR,0,buffer)!=string::npos;
            // Fixing REF and ALT character at POS, if unknown.
            c=ref_sequences.at(record.chrom).at(record.pos-1);
            if (toupper(record.ref.at(0))==VcfReader::UNKNOWN_BASE) record.ref.at(0)=c;
            alt_length=record.alt.length();
            if (toupper(record.alt.at(0))==VcfReader::UNKNOWN_BASE) record.alt.at(0)=c;
            if (toupper(record.alt.at(alt_length-1))==VcfReader::UNKNOWN_BASE) record.alt.at(alt_length-1)=c;
            if (!record.is_breakend_single() && !has_mate) record.info+=VcfReader::INFO_SEPARATOR+MATEID_STR+"="+NEW_MATE_PREFIX+record.id;
            record.print(outstream);
            outstream << "\n";
            // Writing a new mate record if necessary
            if (!record.is_breakend_single() && !has_mate) {
                pos2=record.get_breakend_pos();
                outstream << chromosome2 << VcfReader::VCF_SEPARATOR << std::to_string(pos2) << VcfReader::VCF_SEPARATOR << NEW_MATE_PREFIX << record.id << VcfReader::VCF_SEPARATOR << ref_sequences.at(chromosome2).at(pos2-1) << VcfReader::VCF_SEPARATOR;
                record.get_alt_of_breakend_mate(ref_sequences,buffer);
                outstream << buffer << VcfReader::VCF_SEPARATOR << (record.qual==-1?VcfReader::VCF_MISSING_STRING_1:std::to_string(record.qual)) << VcfReader::VCF_SEPARATOR << record.filter << VcfReader::VCF_SEPARATOR << VcfReader::SVTYPE_STR << "=" << VcfReader::BND_STR << VcfReader::INFO_SEPARATOR << MATEID_STR << "=" << record.id << VcfReader::VCF_SEPARATOR << record.format;
                for (const auto& item: record.genotypes) outstream << VcfReader::VCF_SEPARATOR << item;
                outstream << "\n";
            }
        }
        else if (record.sv_length<min_sv_length) { record.print(outstream); outstream << "\n"; }
        else {
            chromosome_length=(int32_t)ref_sequences.at(record.chrom).length();
            if (record.sv_type==VcfReader::TYPE_DELETION) {
                if (record.pos+record.sv_length+1<=chromosome_length) {
                    // Left side
                    c=ref_sequences.at(record.chrom).at(record.pos-1);
                    outstream << record.chrom << VcfReader::VCF_SEPARATOR << std::to_string(record.pos) << VcfReader::VCF_SEPARATOR << record.id << VcfReader::VCF_SEPARATOR;
                    outstream << c << VcfReader::VCF_SEPARATOR << c << VcfReader::BREAKEND_CHAR_OPEN << record.chrom << VcfReader::BREAKEND_SEPARATOR << std::to_string(record.pos+record.sv_length+1) << VcfReader::BREAKEND_CHAR_OPEN << VcfReader::VCF_SEPARATOR;
                    outstream << (record.qual==-1?VcfReader::VCF_MISSING_STRING_1:std::to_string(record.qual)) << VcfReader::VCF_SEPARATOR << record.filter << VcfReader::VCF_SEPARATOR;
                    outstream << VcfReader::SVTYPE_STR << "=" << VcfReader::BND_STR << VcfReader::INFO_SEPARATOR << MATEID_STR << "=" << NEW_MATE_PREFIX << record.id << VcfReader::VCF_SEPARATOR;
                    outstream << record.format;
                    for (const auto& item: record.genotypes) outstream << VcfReader::VCF_SEPARATOR << item;
                    outstream << "\n";
                    // Right side
                    c=ref_sequences.at(record.chrom).at(record.pos+record.sv_length+1-1);
                    outstream << record.chrom << VcfReader::VCF_SEPARATOR << std::to_string(record.pos+record.sv_length+1) << VcfReader::VCF_SEPARATOR;
                    outstream << NEW_MATE_PREFIX << record.id << VcfReader::VCF_SEPARATOR;
                    outstream << c << VcfReader::VCF_SEPARATOR;
                    outstream << VcfReader::BREAKEND_CHAR_CLOSE << record.chrom << VcfReader::BREAKEND_SEPARATOR << std::to_string(record.pos) << VcfReader::BREAKEND_CHAR_CLOSE << c << VcfReader::VCF_SEPARATOR;
                    outstream << (record.qual==-1?VcfReader::VCF_MISSING_STRING_1:std::to_string(record.qual)) << VcfReader::VCF_SEPARATOR << record.filter << VcfReader::VCF_SEPARATOR;
                    outstream << VcfReader::SVTYPE_STR << "=" << VcfReader::BND_STR << VcfReader::INFO_SEPARATOR << MATEID_STR << "=" << record.id << VcfReader::VCF_SEPARATOR;
                    outstream << record.format;
                    for (const auto& item: record.genotypes) outstream << VcfReader::VCF_SEPARATOR << item;
                    outstream << "\n";
                }
                else { /* Filtered out: this would be a single breakend with no inserted sequence. */ }
            }
            else if (record.sv_type==VcfReader::TYPE_INVERSION) {
                if (record.pos+record.sv_length>chromosome_length) record.sv_length=(int32_t)(chromosome_length-(record.pos+1));
                // Adjacency 1, left side.
                c=ref_sequences.at(record.chrom).at(record.pos-1);
                outstream << record.chrom << VcfReader::VCF_SEPARATOR << std::to_string(record.pos) << VcfReader::VCF_SEPARATOR << record.id << VcfReader::VCF_SEPARATOR;
                outstream << c << VcfReader::VCF_SEPARATOR << c << VcfReader::BREAKEND_CHAR_CLOSE << record.chrom << VcfReader::BREAKEND_SEPARATOR << std::to_string(record.pos+record.sv_length) << VcfReader::BREAKEND_CHAR_CLOSE << VcfReader::VCF_SEPARATOR;
                outstream << (record.qual==-1?VcfReader::VCF_MISSING_STRING_1:std::to_string(record.qual)) << VcfReader::VCF_SEPARATOR << record.filter << VcfReader::VCF_SEPARATOR;
                outstream << VcfReader::SVTYPE_STR << "=" << VcfReader::BND_STR << VcfReader::INFO_SEPARATOR << MATEID_STR << "=" << NEW_MATE_PREFIX << record.id << VcfReader::VCF_SEPARATOR;
                outstream << record.format;
                for (const auto& item: record.genotypes) outstream << VcfReader::VCF_SEPARATOR << item;
                outstream << "\n";
                // Adjacency 1, right side.
                c=ref_sequences.at(record.chrom).at(record.pos+record.sv_length-1);
                outstream << record.chrom << VcfReader::VCF_SEPARATOR;
                outstream << std::to_string(record.pos+record.sv_length) << VcfReader::VCF_SEPARATOR << NEW_MATE_PREFIX << record.id << VcfReader::VCF_SEPARATOR;
                outstream << c << VcfReader::VCF_SEPARATOR << c << VcfReader::BREAKEND_CHAR_CLOSE << record.chrom << VcfReader::BREAKEND_SEPARATOR << std::to_string(record.pos) << VcfReader::BREAKEND_CHAR_CLOSE << VcfReader::VCF_SEPARATOR;
                outstream << (record.qual==-1?VcfReader::VCF_MISSING_STRING_1:std::to_string(record.qual)) << VcfReader::VCF_SEPARATOR << record.filter << VcfReader::VCF_SEPARATOR;
                outstream << VcfReader::SVTYPE_STR << "=" << VcfReader::BND_STR << VcfReader::INFO_SEPARATOR << MATEID_STR << "=" << record.id << VcfReader::VCF_SEPARATOR;
                outstream << record.format;
                for (const auto& item: record.genotypes) outstream << VcfReader::VCF_SEPARATOR << item;
                outstream << "\n";
                if (record.pos+record.sv_length+1<=chromosome_length) {
                    // Adjacency 2, left side.
                    c=ref_sequences.at(record.chrom).at(record.pos+1-1);
                    outstream << record.chrom << VcfReader::VCF_SEPARATOR;
                    outstream << std::to_string(record.pos+1) << VcfReader::VCF_SEPARATOR;
                    outstream << record.id << SECOND_ADJACENCY_SUFFIX << VcfReader::VCF_SEPARATOR;
                    outstream << c << VcfReader::VCF_SEPARATOR;
                    outstream << VcfReader::BREAKEND_CHAR_OPEN << record.chrom << VcfReader::BREAKEND_SEPARATOR << std::to_string(record.pos+record.sv_length+1) << VcfReader::BREAKEND_CHAR_OPEN << c << VcfReader::VCF_SEPARATOR;
                    outstream << (record.qual==-1?VcfReader::VCF_MISSING_STRING_1:std::to_string(record.qual)) << VcfReader::VCF_SEPARATOR << record.filter << VcfReader::VCF_SEPARATOR;
                    outstream << VcfReader::SVTYPE_STR << "=" << VcfReader::BND_STR << VcfReader::INFO_SEPARATOR << MATEID_STR << "=" << NEW_MATE_PREFIX << record.id << SECOND_ADJACENCY_SUFFIX << VcfReader::VCF_SEPARATOR;
                    outstream << record.format;
                    for (const auto& item: record.genotypes) outstream << VcfReader::VCF_SEPARATOR << item;
                    outstream << "\n";
                    // Adjacency 2, right side.
                    c=ref_sequences.at(record.chrom).at(record.pos+record.sv_length+1-1);
                    outstream << record.chrom << VcfReader::VCF_SEPARATOR;
                    outstream << std::to_string(record.pos+record.sv_length+1) << VcfReader::VCF_SEPARATOR;
                    outstream << NEW_MATE_PREFIX << record.id << SECOND_ADJACENCY_SUFFIX << VcfReader::VCF_SEPARATOR;
                    outstream << c << VcfReader::VCF_SEPARATOR;
                    outstream << VcfReader::BREAKEND_CHAR_OPEN << record.chrom << VcfReader::BREAKEND_SEPARATOR << std::to_string(record.pos+1) << VcfReader::BREAKEND_CHAR_OPEN << c << VcfReader::VCF_SEPARATOR;
                    outstream << (record.qual==-1?VcfReader::VCF_MISSING_STRING_1:std::to_string(record.qual)) << VcfReader::VCF_SEPARATOR << record.filter << VcfReader::VCF_SEPARATOR;
                    outstream << VcfReader::SVTYPE_STR << "=" << VcfReader::BND_STR << VcfReader::INFO_SEPARATOR << MATEID_STR << "=" << record.id << SECOND_ADJACENCY_SUFFIX << VcfReader::VCF_SEPARATOR;
                    outstream << record.format;
                    for (const auto& item: record.genotypes) outstream << VcfReader::VCF_SEPARATOR << item;
                    outstream << "\n";
                }
            }
            else if (record.sv_type==VcfReader::TYPE_DUPLICATION) {
                if (record.pos+record.sv_length>chromosome_length) record.sv_length=(int32_t)(chromosome_length-(record.pos+1));
                // Left side
                c=ref_sequences.at(record.chrom).at(record.pos+1-1);
                outstream << record.chrom << VcfReader::VCF_SEPARATOR;
                outstream << std::to_string(record.pos+1) << VcfReader::VCF_SEPARATOR;
                outstream << record.id << VcfReader::VCF_SEPARATOR;
                outstream << c << VcfReader::VCF_SEPARATOR;
                outstream << VcfReader::BREAKEND_CHAR_CLOSE << record.chrom << VcfReader::BREAKEND_SEPARATOR << std::to_string(record.pos+record.sv_length) << VcfReader::BREAKEND_CHAR_CLOSE << c << VcfReader::VCF_SEPARATOR;
                outstream << (record.qual==-1?VcfReader::VCF_MISSING_STRING_1:std::to_string(record.qual)) << VcfReader::VCF_SEPARATOR << record.filter << VcfReader::VCF_SEPARATOR;
                outstream << VcfReader::SVTYPE_STR << "=" << VcfReader::BND_STR << VcfReader::INFO_SEPARATOR << MATEID_STR << "=" << NEW_MATE_PREFIX << record.id << VcfReader::VCF_SEPARATOR;
                outstream << record.format;
                for (const auto& item: record.genotypes) outstream << VcfReader::VCF_SEPARATOR << item;
                outstream << "\n";
                // Right side
                c=ref_sequences.at(record.chrom).at(record.pos+record.sv_length-1);
                outstream << record.chrom << VcfReader::VCF_SEPARATOR;
                outstream << std::to_string(record.pos+record.sv_length) << VcfReader::VCF_SEPARATOR;
                outstream << NEW_MATE_PREFIX << record.id << VcfReader::VCF_SEPARATOR;
                outstream << c << VcfReader::VCF_SEPARATOR;
                outstream << c << VcfReader::BREAKEND_CHAR_OPEN << record.chrom << VcfReader::BREAKEND_SEPARATOR << std::to_string(record.pos+1) << VcfReader::BREAKEND_CHAR_OPEN << VcfReader::VCF_SEPARATOR;
                outstream << (record.qual==-1?VcfReader::VCF_MISSING_STRING_1:std::to_string(record.qual)) << VcfReader::VCF_SEPARATOR << record.filter << VcfReader::VCF_SEPARATOR;
                outstream << VcfReader::SVTYPE_STR << "=" << VcfReader::BND_STR << VcfReader::INFO_SEPARATOR << MATEID_STR << "=" << record.id << VcfReader::VCF_SEPARATOR;
                outstream << record.format;
                for (const auto& item: record.genotypes) outstream << VcfReader::VCF_SEPARATOR << item;
                outstream << "\n";
            }
            else if (record.sv_type==VcfReader::TYPE_INSERTION && !record.is_symbolic) {
                // Adjacency 1
                if (record.pos>0) {
                    c=ref_sequences.at(record.chrom).at(record.pos-1);
                    outstream << record.chrom << VcfReader::VCF_SEPARATOR << std::to_string(record.pos) << VcfReader::VCF_SEPARATOR << record.id << VcfReader::VCF_SEPARATOR;
                    outstream << c << VcfReader::VCF_SEPARATOR;
                    outstream << c << record.alt.substr(1,single_breakend_length) << VcfReader::VCF_MISSING_CHAR_1 << VcfReader::VCF_SEPARATOR;
                    outstream << (record.qual==-1?VcfReader::VCF_MISSING_STRING_1:std::to_string(record.qual)) << VcfReader::VCF_SEPARATOR << record.filter << VcfReader::VCF_SEPARATOR;
                    outstream << VcfReader::SVTYPE_STR << "=" << VcfReader::BND_STR << VcfReader::INFO_SEPARATOR << VcfReader::VCF_SEPARATOR;
                    outstream << record.format;
                    for (const auto& item: record.genotypes) outstream << VcfReader::VCF_SEPARATOR << item;
                    outstream << "\n";
                }
                // Adjacency 2
                if (record.pos+1<=chromosome_length) {
                    c=ref_sequences.at(record.chrom).at(record.pos+1-1);
                    outstream << record.chrom << VcfReader::VCF_SEPARATOR;
                    outstream << std::to_string(record.pos+1) << VcfReader::VCF_SEPARATOR;
                    outstream << NEW_MATE_PREFIX << record.id << VcfReader::VCF_SEPARATOR;
                    outstream << c << VcfReader::VCF_SEPARATOR;
                    outstream << VcfReader::VCF_MISSING_CHAR_1 << record.alt.substr(record.alt.length()-single_breakend_length,single_breakend_length) << c << VcfReader::VCF_SEPARATOR;
                    outstream << (record.qual==-1?VcfReader::VCF_MISSING_STRING_1:std::to_string(record.qual)) << VcfReader::VCF_SEPARATOR << record.filter << VcfReader::VCF_SEPARATOR;
                    outstream << VcfReader::SVTYPE_STR << "=" << VcfReader::BND_STR << VcfReader::INFO_SEPARATOR << VcfReader::VCF_SEPARATOR;
                    outstream << record.format;
                    for (const auto& item: record.genotypes) outstream << VcfReader::VCF_SEPARATOR << item;
                    outstream << "\n";
                }
            }
            else { record.print(outstream); outstream << "\n"; }
        }
    });
    outstream.close();
    cerr << "done.\n";
    return 0;
}
