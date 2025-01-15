#include "misc.hpp"
#include "VcfReader.hpp"

#include <cstdlib>
#include <iostream>

using std::stoul;
using std::cerr;
using std::cout;
using std::to_string;
using std::ofstream;
using sv_merge::run_command;
using sv_merge::VcfRecord;
using sv_merge::VcfReader;


ofstream outstream;
unordered_set<int32_t> tmp_set;
pair<double,double> tmp_pair;
string tmp_string;

void test_callback(VcfRecord& record) { record.print(outstream); outstream << '\n'; }

// TODO: implement proper test
void test_callback_2(VcfRecord& record) {
    record.get_samples_with_alt(tmp_set);

    vector<int32_t> tmp_vector;
    tmp_vector.reserve(tmp_set.size());
    for (auto i: tmp_set) tmp_vector.push_back(i);
    sort(tmp_vector.begin(),tmp_vector.end());
    for (size_t i=0; i<tmp_vector.size(); i++) { outstream << to_string(tmp_vector.at(i)) << ","; }
    outstream << "\n";
}

/**
 * For PBSV.
 * Remark: SIGMA_MULTIPLE must be set to 1 in get_confidence_interval() for this test to work.
 */
void test_callback_3(VcfRecord& record) {
    if (record.get_info_field(VcfReader::CIPOS_STR,0,tmp_string)==string::npos) outstream << ".,\n";
    else {
        record.get_confidence_interval_pos(tmp_pair);
        outstream << tmp_pair.first << "," << tmp_pair.second << ",\n";
    }
}

/**
 * For Sniffles.
 * Remark: SIGMA_MULTIPLE must be set to 1 in get_confidence_interval() for this test to work.
 */
void test_callback_4(VcfRecord& record) {
    if (record.get_info_field(VcfReader::STDEV_POS_STR,0,tmp_string)==string::npos) outstream << ".,";
    else {
        record.get_confidence_interval_pos(tmp_pair);
        outstream << tmp_pair.second << ",";
    }
    if (record.get_info_field(VcfReader::STDEV_LEN_STR,0,tmp_string)==string::npos) outstream << ".,\n";
    else {
        record.get_confidence_interval_length(tmp_pair);
        outstream << tmp_pair.second << ",\n";
    }
}


/**
 * Performs `input_vcf -> test_vcf` and compares `test_vcf` to `truth_vcf`.
 */
void test(const function<void(VcfRecord& record)>& callback, const string& test_id, const string& chromosome, bool high_qual_only, double min_qual, bool pass_only, int32_t min_sv_length, int32_t n_samples_to_load, double min_af, double min_nmf, const path& input_vcf, const path& truth_vcf, const path& test_vcf) {
    cerr << "Test ID: " << test_id << '\n';
    cerr << "   chromosome: " << chromosome << '\n';
    cerr << "   high_qual_only: " << high_qual_only << '\n';
    cerr << "   min_qual: " << min_qual << '\n';
    cerr << "   pass_only: " << pass_only << '\n';
    cerr << "   min_sv_length: " << min_sv_length << '\n';
    cerr << "   n_samples_to_load: " << n_samples_to_load << '\n';
    cerr << "   min_allele_frequency: " << min_af << '\n';
    cerr << "   min_nonmissing_frequency: " << min_nmf << '\n';
    VcfReader reader(input_vcf,10000,high_qual_only,min_qual,pass_only,min_sv_length,INT32_MAX,n_samples_to_load,min_af,min_nmf);
    outstream.open(test_vcf);
    reader.for_record_in_vcf(callback);
    outstream.close();
    for (const auto& sample_id: reader.sample_ids) cerr << sample_id << "  ";
    cerr << "\n";
    string command = "diff --brief "+test_vcf.string()+" "+truth_vcf.string();
    run_command(command);
}


void test_single_sample_vcf(const path& caller_vcf, const bool filter_by_qual, const string& caller_id, const vector<string>& CHROMOSOMES, const vector<double>& MIN_QUALS, const vector<size_t>& MIN_SV_LENGTHS, const path& input_vcf, const path& truth_vcf, const path& test_vcf) {
    bool high_qual_only, pass_only;
    size_t min_sv_length;
    double min_qual, min_af, min_nmf;
    string command;

    for (auto& chromosome: CHROMOSOMES) {
        command=R"(bcftools filter --include "ALT!=\"<INS>\" && ALT!=\"<CNV>\"" --regions )"+chromosome+" "+caller_vcf.string();
        run_command(command,input_vcf);
        // Testing get_confidence_interval_*()
        high_qual_only=false; pass_only=false; min_sv_length=0; min_af=0.0; min_nmf=0.0;
        if (caller_id=="PBSV") {
            command=R"(bcftools query -f '%CIPOS,\n' )"+input_vcf.string();
            run_command(command,truth_vcf);
            test(test_callback_3,caller_id,chromosome,high_qual_only,min_qual,pass_only,min_sv_length,1,min_af,min_nmf,input_vcf,truth_vcf,test_vcf);
        }
        else if (caller_id=="SNIFFLES") {
            command=R"(bcftools query -f '%STDEV_POS,%STDEV_LEN,\n' )"+input_vcf.string();
            run_command(command,truth_vcf);
            test(test_callback_4,caller_id,chromosome,high_qual_only,min_qual,pass_only,min_sv_length,1,min_af,min_nmf,input_vcf,truth_vcf,test_vcf);
        }
        else if (caller_id=="PAV") {
            // NOP: PAV has no confidence interval information.
        }
        // Filtering by QUAL
        if (filter_by_qual) {
            high_qual_only=true; pass_only=false; min_sv_length=0; min_af=0.0; min_nmf=0.0;
            for (auto& min_qual: MIN_QUALS) {
                command=R"(bcftools filter --include "QUAL>=)"+to_string(min_qual)+R"(" )"+input_vcf.string()+R"( | grep "^[^#]")";
                run_command(command,truth_vcf);
                test(test_callback,caller_id,chromosome,high_qual_only,min_qual,pass_only,min_sv_length,1,min_af,min_nmf,input_vcf,truth_vcf,test_vcf);
            }
        }
        // Filtering by FILTER
        high_qual_only=false; pass_only=true; min_sv_length=0; min_af=0.0; min_nmf=0.0;
        command=R"(bcftools filter --include "FILTER=\"PASS\"" )"+input_vcf.string()+R"( | grep "^[^#]")";
        run_command(command,truth_vcf);
        min_qual=0.0;
        test(test_callback,caller_id,chromosome,high_qual_only,min_qual,pass_only,min_sv_length,1,min_af,min_nmf,input_vcf,truth_vcf,test_vcf);
        // Filtering by SVLEN
        for (auto& min_sv_length: MIN_SV_LENGTHS) {
            high_qual_only=false; pass_only=false; min_af=0.0; min_nmf=0.0;
            command=R"(bcftools filter --include "SVLEN>=)"+to_string(min_sv_length)+" || SVLEN<=-"+to_string(min_sv_length)+R"( || SVTYPE==\"BND\"" )"+input_vcf.string()+R"( | grep "^[^#]")";
            run_command(command,truth_vcf);
            test(test_callback,caller_id,chromosome,high_qual_only,min_qual,pass_only,min_sv_length,1,min_af,min_nmf,input_vcf,truth_vcf,test_vcf);
        }
        // Removing sample columns
        high_qual_only=false; min_qual=0.0; pass_only=false; min_sv_length=0; min_af=0.0; min_nmf=0.0;
        command="cut -f 1-9 "+input_vcf.string()+R"( | grep "^[^#]")";
        run_command(command,truth_vcf);
        test(test_callback,caller_id,chromosome,high_qual_only,min_qual,pass_only,min_sv_length,0,min_af,min_nmf,input_vcf,truth_vcf,test_vcf);
    }
}


/**
 * Currently tuned for the "raw" HPRC VCF.
 */
void test_joint_vcf_hprc(const path& joint_vcf, const vector<string>& CHROMOSOMES, const vector<double>& MIN_QUALS, const vector<double>& MIN_AFS, const vector<double>& MIN_NMFS, const string& N_THREADS, const path& tmp1_vcf, const path& tmp2_vcf, const path& input_vcf, const path& truth_vcf, const path& test_vcf) {
    bool high_qual_only, pass_only;
    size_t min_sv_length, n_samples_to_load, n_haplotypes;
    double min_qual, min_af, min_nmf;
    string command;

    n_samples_to_load=45;
    for (auto& chromosome: CHROMOSOMES) {
        // Removing the first sample since it's a haploid reference
        command="bcftools view --threads "+N_THREADS+" "+joint_vcf.string()+" "+chromosome+" | cut -f 1-9,11-$((9+"+to_string(n_samples_to_load)+"))";
        run_command(command,tmp1_vcf);
        n_samples_to_load--; n_haplotypes=n_samples_to_load<<1;
        command="bcftools norm --threads "+N_THREADS+" --multiallelics - "+tmp1_vcf.string();
        run_command(command,tmp2_vcf);
        command="rm -f tmp1_vcf"; run_command(command);
        command="bcftools annotate --threads "+N_THREADS+" -x INFO/AC,INFO/AF,INFO/AN "+tmp2_vcf.string();
        run_command(command,tmp1_vcf);
        command="rm -f "+tmp2_vcf.string(); run_command(command);
        command="bcftools +fill-tags "+tmp1_vcf.string()+" -Ob -o "+input_vcf.string()+" -- -t AC,AF,AN,F_MISSING"; run_command(command);
        command="rm -f tmp1_vcf"; run_command(command);
        // Testing get_samples_with_alt()
        command="bcftools view -H "+input_vcf.string()+R"( | awk '{for (i=10; i<=54; i++) { if ($i=="0/1" || $i=="0|1" || $i=="1/0" || $i=="1|0" || $i=="1/1" || $i=="1|1" || $i=="./1" || $i==".|1" || $i=="1/." || $i=="1|.") printf("%d,",(i-10)); } printf("\n");}')";
        run_command(command,truth_vcf);
        high_qual_only=false; min_qual=0; pass_only=false; min_sv_length=0; min_af=0.0; min_nmf=0.0;
        test(test_callback_2,"HPRC",chromosome,high_qual_only,min_qual,pass_only,min_sv_length,n_samples_to_load,min_af,min_nmf,input_vcf,truth_vcf,test_vcf);
        // Filtering by QUAL
        high_qual_only=true; pass_only=false; min_sv_length=0; min_af=0.0; min_nmf=0.0;
        for (auto& min_qual: MIN_QUALS) {
            command="bcftools filter --threads "+N_THREADS+R"( --include "QUAL>=)"+to_string(min_qual)+R"(" )"+input_vcf.string()+R"( | grep "^[^#]" | sed 's/PASS/./g')";
            run_command(command,truth_vcf);
            test(test_callback,"HPRC",chromosome,high_qual_only,min_qual,pass_only,min_sv_length,n_samples_to_load,min_af,min_nmf,input_vcf,truth_vcf,test_vcf);
        }
        // Filtering by FILTER skipped, since FILTER is always '.'
        // Filtering by SVLEN skipped, since the SVLEN tag is missing, bcftools does not work for adding it, and
        // `truvari anno svinfo` can add it but alters the VCF in other ways that break the test.
        // Removing all sample columns
        high_qual_only=false; min_qual=0.0; pass_only=false; min_sv_length=0; min_af=0.0; min_nmf=0.0;
        command="cut -f 1-9 "+input_vcf.string()+R"( | grep "^[^#]")";
        run_command(command,truth_vcf);
        test(test_callback,"HPRC",chromosome,high_qual_only,min_qual,pass_only,min_sv_length,0,min_af,min_nmf,input_vcf,truth_vcf,test_vcf);
        // Filtering by AF
        high_qual_only=false; min_qual=0.0; pass_only=false; min_sv_length=0; min_nmf=0.0;
        for (auto& min_af: MIN_AFS) {
            command="bcftools filter --threads "+N_THREADS+R"( --include "AC/)"+to_string(n_haplotypes)+">="+to_string(min_af)+R"(" )"+input_vcf.string()+R"( | grep "^[^#]" | sed 's/PASS/./g')";
            run_command(command,truth_vcf);
            test(test_callback,"HPRC",chromosome,high_qual_only,min_qual,pass_only,min_sv_length,n_samples_to_load,min_af,min_nmf,input_vcf,truth_vcf,test_vcf);
        }
        // Filtering by MIN_NMF
        high_qual_only=false; min_qual=0.0; pass_only=false; min_sv_length=0; min_af=0.0;
        for (auto& min_nmf: MIN_NMFS) {
            command="bcftools filter --threads "+N_THREADS+R"( --include "AN/)"+to_string(n_haplotypes)+">="+to_string(min_nmf)+R"(" )"+input_vcf.string()+R"( | grep "^[^#]" | sed 's/PASS/./g')";
            run_command(command,truth_vcf);
            test(test_callback,"HPRC",chromosome,high_qual_only,min_qual,pass_only,min_sv_length,n_samples_to_load,min_af,min_nmf,input_vcf,truth_vcf,test_vcf);
        }
    }
}


/**
 * Currently tuned for the conventions of a specific file
 */
void test_joint_snp_vcf_hprc(const path& joint_vcf, const vector<string>& CHROMOSOMES, const vector<double>& MIN_QUALS, const vector<double>& MIN_AFS, const vector<double>& MIN_NMFS, const string& N_THREADS, const path& tmp1_vcf, const path& tmp2_vcf, const path& input_vcf, const path& truth_vcf, const path& test_vcf) {
    bool high_qual_only, pass_only;
    size_t min_sv_length, n_samples_to_load, n_haplotypes;
    double min_qual, min_af, min_nmf;
    string command;

    n_samples_to_load=47;
    for (auto& chromosome: CHROMOSOMES) {
        n_haplotypes=n_samples_to_load<<1;
        command="bcftools norm --threads "+N_THREADS+" --multiallelics - "+joint_vcf.string();
        run_command(command,tmp2_vcf);
        command="bcftools annotate --threads "+N_THREADS+" -x INFO/AC,INFO/AF,INFO/AN "+tmp2_vcf.string();
        run_command(command,tmp1_vcf);
        command="rm -f "+tmp2_vcf.string(); run_command(command);
        command="bcftools +fill-tags "+tmp1_vcf.string()+" -Ob -o "+input_vcf.string()+" -- -t AC,AF,AN,F_MISSING"; run_command(command);
        command="rm -f tmp1_vcf"; run_command(command);
        // Testing get_samples_with_alt()
        command="bcftools view -H "+input_vcf.string()+R"( | awk '{for (i=10; i<=56; i++) { if (substr($i,1,3)=="0/1" || substr($i,1,3)=="0|1" || substr($i,1,3)=="1/0" || substr($i,1,3)=="1|0" || substr($i,1,3)=="1/1" || substr($i,1,3)=="1|1" || substr($i,1,3)=="./1" || substr($i,1,3)==".|1" || substr($i,1,3)=="1/." || substr($i,1,3)=="1|.") printf("%d,",(i-10)); } printf("\n");}')";
        run_command(command,truth_vcf);
        high_qual_only=false; min_qual=0; pass_only=false; min_sv_length=0; min_af=0.0; min_nmf=0.0;
        test(test_callback_2,"HPRC",chromosome,high_qual_only,min_qual,pass_only,min_sv_length,n_samples_to_load,min_af,min_nmf,input_vcf,truth_vcf,test_vcf);
        // Filtering by QUAL
        high_qual_only=true; pass_only=false; min_sv_length=0; min_af=0.0; min_nmf=0.0;
        for (auto& min_qual: MIN_QUALS) {
            command="bcftools filter --threads "+N_THREADS+R"( --include "QUAL>=)"+to_string(min_qual)+R"(" )"+input_vcf.string()+R"( | grep "^[^#]" | sed 's/PASS/./g')";
            run_command(command,truth_vcf);
            test(test_callback,"HPRC",chromosome,high_qual_only,min_qual,pass_only,min_sv_length,n_samples_to_load,min_af,min_nmf,input_vcf,truth_vcf,test_vcf);
        }
        // Filtering by FILTER skipped, since FILTER is always '.'
        // Removing all sample columns
        high_qual_only=false; min_qual=0.0; pass_only=false; min_sv_length=0; min_af=0.0; min_nmf=0.0;
        command="cut -f 1-9 "+input_vcf.string()+R"( | grep "^[^#]")";
        run_command(command,truth_vcf);
        test(test_callback,"HPRC",chromosome,high_qual_only,min_qual,pass_only,min_sv_length,0,min_af,min_nmf,input_vcf,truth_vcf,test_vcf);
        // Filtering by AF
        high_qual_only=false; min_qual=0.0; pass_only=false; min_sv_length=0; min_nmf=0.0;
        for (auto& min_af: MIN_AFS) {
            command="bcftools filter --threads "+N_THREADS+R"( --include "AC/)"+to_string(n_haplotypes)+">="+to_string(min_af)+R"(" )"+input_vcf.string()+R"( | grep "^[^#]" | sed 's/PASS/./g')";
            run_command(command,truth_vcf);
            test(test_callback,"HPRC",chromosome,high_qual_only,min_qual,pass_only,min_sv_length,n_samples_to_load,min_af,min_nmf,input_vcf,truth_vcf,test_vcf);
        }
        // Filtering by MIN_NMF
        high_qual_only=false; min_qual=0.0; pass_only=false; min_sv_length=0; min_af=0.0;
        for (auto& min_nmf: MIN_NMFS) {
            command="bcftools filter --threads "+N_THREADS+R"( --include "AN/)"+to_string(n_haplotypes)+">="+to_string(min_nmf)+R"(" )"+input_vcf.string()+R"( | grep "^[^#]" | sed 's/PASS/./g')";
            run_command(command,truth_vcf);
            test(test_callback,"HPRC",chromosome,high_qual_only,min_qual,pass_only,min_sv_length,n_samples_to_load,min_af,min_nmf,input_vcf,truth_vcf,test_vcf);
        }
    }
}


int main(int argc, char* argv[]) {
    const path ROOT_DIR = path(argv[1]);  // Assumed to contain every raw input file used for testing
    const char* BCFTOOLS_PLUGINS_DIR = argv[2];  // Needed by bcftools +fill-tags
    const string N_THREADS = string(argv[3]);

    /**
     * Testing ranges
     */
    const vector<size_t> MIN_SV_LENGTHS = {0, 50, 100, 500, 1000};
    const vector<double> MIN_QUALS = {0.0, 10.0, 20.0, 30.0};
    const vector<double> MIN_AFS = {0.02, 0.08, 0.16, 0.32, 0.64};
    const vector<double> MIN_NMFS = {0.02, 0.08, 0.16, 0.32, 0.64};
    // sniffles (and maybe others) doesn't have chrM calls.
    const vector<string> CHROMOSOMES_SINGLE_SAMPLE = {"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"};
    // chrY and chrM are not present in the HPRC VCF. Our AF filters in chrX behave differently from bcftools and would
    // require a better implementation of `test_callback()`.
    const vector<string> CHROMOSOMES_HPRC = {"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"};
    const vector<string> SAMPLE_IDS = {"NA24385", "HG03125", "HG00512"};
    const string HPRC_FILE_ID = "hprc-v1.1-mc-chm13.raw.sv";
    const string HPRC_SNP_FILE_ID = "HPRC.DV.joint.g";

    string command;
    path input_vcf, test_vcf, truth_vcf, tmp1_vcf, tmp2_vcf;  // Temporary files created during the tests

    setenv("BCFTOOLS_PLUGINS",BCFTOOLS_PLUGINS_DIR,1);
    input_vcf=ROOT_DIR/"input.vcf";
    test_vcf=ROOT_DIR/"test.vcf";
    truth_vcf=ROOT_DIR/"truth.vcf";
    tmp1_vcf=ROOT_DIR/"tmp1.vcf";
    tmp2_vcf=ROOT_DIR/"tmp2.vcf";
    for (auto& sample_id: SAMPLE_IDS) {
        test_single_sample_vcf(ROOT_DIR/(sample_id+".sniffles.vcf.gz"),true,"SNIFFLES",CHROMOSOMES_SINGLE_SAMPLE,MIN_QUALS,MIN_SV_LENGTHS,input_vcf,truth_vcf,test_vcf);
        test_single_sample_vcf(ROOT_DIR/(sample_id+".pbsv.vcf.gz"),false/*Filtering by QUAL skipped, since QUAL is always '.'*/,"PBSV",CHROMOSOMES_SINGLE_SAMPLE,MIN_QUALS,MIN_SV_LENGTHS,input_vcf,truth_vcf,test_vcf);
        test_single_sample_vcf(ROOT_DIR/(sample_id+".pav_sv.vcf.gz"),false/*Filtering by QUAL skipped, since QUAL is always '.'*/,"PAV",CHROMOSOMES_SINGLE_SAMPLE,MIN_QUALS,MIN_SV_LENGTHS,input_vcf,truth_vcf,test_vcf);
    }
    command="rm -f "+input_vcf.string()+" "+test_vcf.string()+" "+truth_vcf.string(); run_command(command);
    test_joint_vcf_hprc(ROOT_DIR/(HPRC_FILE_ID+".vcf.gz"),CHROMOSOMES_HPRC,MIN_QUALS,MIN_AFS,MIN_NMFS,N_THREADS,tmp1_vcf,tmp2_vcf,input_vcf,truth_vcf,test_vcf);
    command="rm -f "+input_vcf.string()+" "+test_vcf.string()+" "+truth_vcf.string()+" "+tmp1_vcf.string()+" "+tmp2_vcf.string(); run_command(command);
    test_joint_snp_vcf_hprc(ROOT_DIR/(HPRC_SNP_FILE_ID+".vcf.gz"),CHROMOSOMES_HPRC,MIN_QUALS,MIN_AFS,MIN_NMFS,N_THREADS,tmp1_vcf,tmp2_vcf,input_vcf,truth_vcf,test_vcf);
    command="rm -f "+input_vcf.string()+" "+test_vcf.string()+" "+truth_vcf.string()+" "+tmp1_vcf.string()+" "+tmp2_vcf.string(); run_command(command);

    return 0;
}