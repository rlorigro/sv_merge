#include "misc.hpp"
#include "VcfReader.hpp"

#include <iostream>

using std::stoul;
using std::cerr;
using std::cout;
using std::to_string;
using std::ofstream;
using sv_merge::run_command;
using sv_merge::VcfRecord;
using sv_merge::VcfReader;


/**
 * Testing ranges
 */
const vector<int> MIN_SV_LENGTHS = {0, 50, 100, 500, 1000};
const vector<float> MIN_QUALS = {0.0, 10.0, 20.0, 30.0};
const vector<float> MIN_AFS = {0.02, 0.08, 0.16, 0.32, 0.64};
const vector<float> MIN_NMFS = {0.02, 0.08, 0.16, 0.32, 0.64};
const vector<string> CHROMOSOMES = {"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"};
const vector<string> SAMPLE_IDS = {"NA24385", "HG03125", "HG00512"};
const string HPRC_FILE_ID = "hprc-v1.1-mc-chm13.raw.sv";
const uint32_t PROGRESS_EVERY_LINES = 10000;  // Arbitrary

string N_THREADS;

/**
 * Assumed to contain every raw input file used for testing
 */
path ROOT_DIR;

/**
 * Temporary files created during the tests
 */
path input_vcf, test_vcf, truth_vcf, tmp1_vcf, tmp2_vcf;
ofstream outstream;


void test_callback(VcfRecord& record) { record.print(outstream); outstream << '\n'; }


/**
 * Performs `input_vcf -> test_vcf` and compares `test_vcf` to `truth_vcf`.
 */
void test(const string& test_id, const string& chromosome, const bool& is_diploid, const bool& high_qual_only, const float& min_qual, const bool& pass_only, const uint32_t min_sv_length, const bool& skip_samples, const uint32_t n_samples_in_vcf, const float& min_af, const float& min_nmf) {
    cerr << "Test ID: " << test_id << '\n';
    cerr << "   chromosome: " << chromosome << '\n';
    cerr << "   is_diploid: " << is_diploid << '\n';
    cerr << "   high_qual_only: " << high_qual_only << '\n';
    cerr << "   min_qual: " << min_qual << '\n';
    cerr << "   pass_only: " << pass_only << '\n';
    cerr << "   min_sv_length: " << min_sv_length << '\n';
    cerr << "   skip_samples: " << skip_samples << '\n';
    cerr << "   n_samples_in_vcf: " << n_samples_in_vcf << '\n';
    cerr << "   min_allele_frequency: " << min_af << '\n';
    cerr << "   min_nonmissing_frequency: " << min_nmf << '\n';
    VcfReader reader(input_vcf,test_callback,PROGRESS_EVERY_LINES,chromosome,is_diploid,high_qual_only,min_qual,pass_only,min_sv_length,skip_samples,n_samples_in_vcf,min_af,min_nmf);
    outstream.open(test_vcf);
    reader.for_record_in_vcf();
    outstream.close();
    string command = "diff --brief "+test_vcf.string()+" "+truth_vcf.string();
    run_command(command);
}


void test_single_sample_vcf(const path& caller_vcf, const bool& filter_by_qual, const string& caller_id) {
    bool is_diploid, high_qual_only, pass_only, skip_samples;
    uint32_t min_sv_length;
    const uint32_t n_samples_in_vcf = 1;
    float min_qual, min_af, min_nmf;
    string command;

    for (auto& chromosome: CHROMOSOMES) {
        is_diploid=chromosome!="chrY"&&chromosome!="chrM";
        command=R"(bcftools filter --include "ALT!=\"<INS>\" && ALT!=\"<CNV>\"" --regions )"+chromosome+" "+caller_vcf.string();
        run_command(command,input_vcf);
        // Filtering by QUAL
        if (filter_by_qual) {
            high_qual_only=true; pass_only=false; min_sv_length=0; skip_samples=false; min_af=0.0; min_nmf=0.0;
            for (auto& min_qual: MIN_QUALS) {
                command=R"(bcftools filter --include "QUAL>=)"+to_string(min_qual)+R"(" )"+input_vcf.string()+R"( | grep "^[^#]")";
                run_command(command,truth_vcf);
                test(caller_id,chromosome,is_diploid,high_qual_only,min_qual,pass_only,min_sv_length,skip_samples,n_samples_in_vcf,min_af,min_nmf);
            }
        }
        // Filtering by FILTER
        high_qual_only=false; pass_only=true; min_sv_length=0; skip_samples=false; min_af=0.0; min_nmf=0.0;
        command=R"(bcftools filter --include "FILTER=\"PASS\"" )"+input_vcf.string()+R"( | grep "^[^#]")";
        run_command(command,truth_vcf);
        min_qual=0.0;
        test(caller_id,chromosome,is_diploid,high_qual_only,min_qual,pass_only,min_sv_length,skip_samples,n_samples_in_vcf,min_af,min_nmf);
        // Filtering by SVLEN
        for (auto& min_sv_length: MIN_SV_LENGTHS) {
            high_qual_only=false; pass_only=false; skip_samples=false; min_af=0.0; min_nmf=0.0;
            command=R"(bcftools filter --include "SVLEN>=)"+to_string(min_sv_length)+" || SVLEN<=-"+to_string(min_sv_length)+R"( || SVTYPE==\"BND\"" )"+input_vcf.string()+R"( | grep "^[^#]")";
            run_command(command,truth_vcf);
            test(caller_id,chromosome,is_diploid,high_qual_only,min_qual,pass_only,min_sv_length,skip_samples,n_samples_in_vcf,min_af,min_nmf);
        }
        // Removing sample columns
        high_qual_only=false; min_qual=0.0; pass_only=false; min_sv_length=0; skip_samples=true; min_af=0.0; min_nmf=0.0;
        command="cut -f 1-9 "+input_vcf.string()+R"( | grep "^[^#]")";
        run_command(command,truth_vcf);
        test(caller_id,chromosome,is_diploid,high_qual_only,min_qual,pass_only,min_sv_length,skip_samples,n_samples_in_vcf,min_af,min_nmf);
    }
}


/**
 * Currently tuned for the "raw" HPRC VCF.
 */
void test_joint_vcf_hprc(const path& joint_vcf) {
    bool is_diploid, high_qual_only, pass_only, skip_samples;
    uint32_t min_sv_length, n_samples_in_vcf, n_haplotypes;
    float min_qual, min_af, min_nmf;
    string command;

    n_samples_in_vcf=45;
    for (auto& chromosome: CHROMOSOMES) {
        is_diploid=chromosome!="chrY"&&chromosome!="chrM";
        // Removing the first sample since it's a haploid reference
        command="bcftools view --threads "+N_THREADS+" "+joint_vcf.string()+" "+chromosome+" | cut -f 1-9,11-$((9+"+to_string(n_samples_in_vcf)+"))";
        run_command(command,tmp1_vcf);
        n_samples_in_vcf--; n_haplotypes=n_samples_in_vcf<<1;
        command="bcftools norm --threads "+N_THREADS+" --multiallelics - "+tmp1_vcf.string();
        run_command(command,tmp2_vcf);
        command="rm -f tmp1_vcf"; run_command(command);
        command="bcftools annotate --threads "+N_THREADS+" -x INFO/AC,INFO/AF,INFO/AN "+tmp2_vcf.string();
        run_command(command,tmp1_vcf);
        command="rm -f "+tmp2_vcf.string(); run_command(command);
        command="bcftools +fill-tags "+tmp1_vcf.string()+" -Ob -o "+input_vcf.string()+" -- -t AC,AF,AN,F_MISSING"; run_command(command);
        command="rm -f tmp1_vcf"; run_command(command);
        // Filtering by QUAL
        high_qual_only=true; pass_only=false; min_sv_length=0; skip_samples=false; min_af=0.0; min_nmf=0.0;
        for (auto& min_qual: MIN_QUALS) {
            command="bcftools filter --threads "+N_THREADS+R"( --include "QUAL>=)"+to_string(min_qual)+R"(" )"+input_vcf.string()+R"( | grep "^[^#]" | sed 's/PASS/./g')";
            run_command(command,truth_vcf);
            test("HPRC",chromosome,is_diploid,high_qual_only,min_qual,pass_only,min_sv_length,skip_samples,n_samples_in_vcf,min_af,min_nmf);
        }
        // Filtering by FILTER skipped, since FILTER is always '.'
        // Filtering by SVLEN skipped, since the SVLEN tag is missing, bcftools does not work for adding it, and
        // `truvari anno svinfo` can add it but alters the VCF in other ways that break the test.
        // Removing all sample columns
        high_qual_only=false; min_qual=0.0; pass_only=false; min_sv_length=0; skip_samples=true; min_af=0.0; min_nmf=0.0;
        command="cut -f 1-9 "+input_vcf.string()+R"( | grep "^[^#]")";
        run_command(command,truth_vcf);
        test("HPRC",chromosome,is_diploid,high_qual_only,min_qual,pass_only,min_sv_length,skip_samples,n_samples_in_vcf,min_af,min_nmf);
        // Filtering by AF
        high_qual_only=false; min_qual=0.0; pass_only=false; min_sv_length=0; skip_samples=false; min_nmf=0.0;
        for (auto& min_af: MIN_AFS) {
            command="bcftools filter --threads "+N_THREADS+R"( --include "AC/)"+to_string(n_haplotypes)+">="+to_string(min_af)+R"(" )"+input_vcf.string()+R"( | grep "^[^#]" | sed 's/PASS/./g')";
            run_command(command,truth_vcf);
            test("HPRC",chromosome,is_diploid,high_qual_only,min_qual,pass_only,min_sv_length,skip_samples,n_samples_in_vcf,min_af,min_nmf);
        }
        // Filtering by MIN_NMF
        high_qual_only=false; min_qual=0.0; pass_only=false; min_sv_length=0; skip_samples=false; min_af=0.0;
        for (auto& min_nmf: MIN_NMFS) {
            command="bcftools filter --threads "+N_THREADS+R"( --include "AN/)"+to_string(n_haplotypes)+">="+to_string(min_nmf)+R"(" )"+input_vcf.string()+R"( | grep "^[^#]" | sed 's/PASS/./g')";
            run_command(command,truth_vcf);
            test("HPRC",chromosome,is_diploid,high_qual_only,min_qual,pass_only,min_sv_length,skip_samples,n_samples_in_vcf,min_af,min_nmf);
        }
    }
}

/**
 * Requirements:
 * - the environment variable `BCFTOOLS_PLUGINS` should be set.
 */
int main(int argc, char* argv[]) {
    ROOT_DIR=path(argv[1]);
    N_THREADS=string(argv[2]);

    string command;

    // Initializing locations of temporary files
    input_vcf=ROOT_DIR/"input.vcf";
    test_vcf=ROOT_DIR/"test.vcf";
    truth_vcf=ROOT_DIR/"truth.vcf";
    tmp1_vcf=ROOT_DIR/"tmp1.vcf";
    tmp2_vcf=ROOT_DIR/"tmp2.vcf";

    for (auto& sample_id: SAMPLE_IDS) {
        test_single_sample_vcf(ROOT_DIR/(sample_id+".sniffles.vcf.gz"),true,"SNIFFLES");
        test_single_sample_vcf(ROOT_DIR/(sample_id+".pbsv.vcf.gz"),false,"PBSV");  // Filtering by QUAL skipped, since QUAL is always '.'
        test_single_sample_vcf(ROOT_DIR/(sample_id+".pav_sv.vcf.gz"),false,"PAV");  // Filtering by QUAL skipped, since QUAL is always '.'
    }
    command="rm -f "+input_vcf.string()+" "+test_vcf.string()+" "+truth_vcf.string(); run_command(command);
    test_joint_vcf_hprc(ROOT_DIR/(HPRC_FILE_ID+".vcf.gz"));
    command="rm -f "+input_vcf.string()+" "+test_vcf.string()+" "+truth_vcf.string()+" "+tmp1_vcf.string()+" "+tmp2_vcf.string(); run_command(command);
    return 0;
}