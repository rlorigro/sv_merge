#include <iostream>
#include <stdexcept>
#include <string>

using std::runtime_error;
using std::string;
using std::cerr;

#include "Filesystem.hpp"

using ghc::filesystem::path;

#include "wfa2/include/wfa2lib/bindings/cpp/WFAligner.hpp"

using namespace wfa;

int main(){
    // Create a WFAligner
    WFAlignerGapAffine aligner(4,6,2,WFAligner::Alignment,WFAligner::MemoryHigh);

    string pattern = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT";
    string text    = "TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
    aligner.alignEnd2End(pattern,text); // Align

    // Display CIGAR & score
    string cigar = aligner.getCIGAR(true);
    cerr << "CIGAR: " << cigar  << '\n';
    cerr << "Alignment score " << aligner.getAlignmentScore() << '\n';
    return 0;
}
