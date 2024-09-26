#include "VcfReader.hpp"
#include "VariantGraph.hpp"

using sv_merge::run_command;
using sv_merge::VcfReader;
using sv_merge::VariantGraph;

using std::max;
using std::min;
using std::to_string;
using std::stoi;

using std::filesystem::path;
using std::filesystem::directory_iterator;
using coord_t = pair<int32_t,int32_t>;


int main(int argc, char* argv[]) {
    const path INPUT_VCF = path(argv[1]);
    const path CHUNKS_DIR = path(argv[2]);
    const int32_t LINES_PER_CHUNK = stoi(argv[3]);
    const string CHR_ID = argv[4];

    const int32_t BUFFER = 100;  // Arbitrary
    size_t i, j;
    int32_t min_pos, max_pos;
    int32_t n_ins, n_del, n_snp, n_replacement;
    int32_t n_reference_nodes, n_nodes, n_edges;
    string current_file, chromosome, command;
    coord_t coords;
    vector<int32_t> first_positions;  // Zero-based
    vector<VcfRecord> records;
    unordered_map<string,string> chromosomes;

    // Creating chunks
    command.clear(); command.append("rm -rf "+CHUNKS_DIR.string()+" ; mkdir "+CHUNKS_DIR.string()); run_command(command);
    command.clear(); command.append("bcftools view --header-only "+INPUT_VCF.string()+" | tail -n 2 > header.txt"); run_command(command);
    command.clear(); command.append("bcftools view --no-header "+INPUT_VCF.string()+" "+CHR_ID+" > "+CHUNKS_DIR.string()+"/tmp.txt"); run_command(command);
    command.clear(); command.append("cd "+CHUNKS_DIR.string()+"; split -l "+to_string(LINES_PER_CHUNK)+" -a 3 tmp.txt; rm -f tmp.txt"); run_command(command);
    command.clear(); command.append("for FILE in $(ls "+CHUNKS_DIR.string()+"/x*); do cat header.txt ${FILE} > ${FILE}.vcf ; done"); run_command(command);

    // Analyzing every chunk
    for (const auto& entry: directory_iterator(CHUNKS_DIR)) {
        current_file.clear(); current_file.append(entry.path().string());
        if (!current_file.ends_with(".vcf")) continue;
        cerr << "Testing " << current_file << "...\n";
        VcfReader reader(current_file);
        records.clear();
        reader.for_record_in_vcf([&](VcfRecord& record) { records.push_back(record); });

        // Removing records with POS=0, for simplicity.
        j=-1;
        for (i=0; i<records.size(); i++) {
            if (records.at(i).pos==0) continue;
            j++;
            VcfRecord tmp_record = records.at(j);
            records.at(j)=records.at(i);
            records.at(i)=tmp_record;
        }
        records.erase(records.begin()+j+1,records.end());

        // Computing statistics, collecting the distinct first position of every reference node, and making all POS
        // values relative.
        first_positions.clear();
        min_pos=INT32_MAX; max_pos=0; n_ins=0; n_del=0; n_snp=0; n_replacement=0;
        for (auto& record: records) {
            min_pos=min(min_pos,record.pos);
            max_pos=max(max_pos,record.pos);
            record.get_reference_coordinates(false,coords);
            if (record.sv_type==VcfReader::TYPE_INSERTION) {
                n_ins++;
                first_positions.emplace_back(coords.first);
            }
            else if (record.sv_type==VcfReader::TYPE_DELETION) {
                n_del++;
                first_positions.emplace_back(coords.first);
                first_positions.emplace_back(coords.second);
            }
            else if (record.sv_type==VcfReader::TYPE_SNP) {
                n_snp++;
                first_positions.emplace_back(coords.first);
                first_positions.emplace_back(coords.second);
            }
            else if (record.sv_type==VcfReader::TYPE_REPLACEMENT) {
                n_replacement++;
                first_positions.emplace_back(coords.first);
                first_positions.emplace_back(coords.second);
            }
        }
        min_pos=max(min_pos-BUFFER,0); max_pos=max_pos+BUFFER;
        for (auto& record: records) record.pos-=min_pos;
        max_pos-=min_pos; min_pos=0;
        std::sort(first_positions.begin(),first_positions.end());
        j=0;
        for (i=1; i<first_positions.size(); i++) {
            if (first_positions.at(i)==first_positions.at(j)) continue;
            j++;
            int32_t tmp = first_positions.at(j);
            first_positions.at(j)=first_positions.at(i);
            first_positions.at(i)=tmp;
        }
        n_reference_nodes=j+2;

        // Building the graph
        chromosome.clear();
        for (i=min_pos; i<=(size_t)max_pos; i++) chromosome.push_back('A');
        chromosomes.clear();
        chromosomes.insert({CHR_ID,chromosome});
        VariantGraph graph(chromosomes);
        graph.build(records,150,15000,false);
        graph.to_gfa("tmp.gfa");
        n_nodes=n_reference_nodes+n_snp+n_replacement+n_ins;
        command.clear(); command.append("grep ^S tmp.gfa | wc -l | xargs > tmp1.txt"); run_command(command);
        command.clear(); command.append("echo "+to_string(n_nodes)+" > tmp2.txt"); run_command(command);
        command.clear(); command.append("diff --brief tmp1.txt tmp2.txt"); run_command(command);
        n_edges=n_reference_nodes-1+n_snp*2+n_replacement*2+n_ins*2+n_del;
        command.clear(); command.append("grep ^L tmp.gfa | wc -l | xargs > tmp1.txt"); run_command(command);
        command.clear(); command.append("echo "+to_string(n_edges)+" > tmp2.txt"); run_command(command);
        command.clear(); command.append("diff --brief tmp1.txt tmp2.txt"); run_command(command);
    }
    return 0;
}