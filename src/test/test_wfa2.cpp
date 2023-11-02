#include <iostream>
#include <stdexcept>
#include <string>

using std::runtime_error;
using std::string;
using std::cerr;

#include "Filesystem.hpp"

using ghc::filesystem::path;

#include "wfa2lib/bindings/cpp/WFAligner.hpp"

using wfa::WFAlignerGapAffine;
using wfa::WFAligner;


int main(){
    // Create a WFAligner
    WFAlignerGapAffine aligner(4,6,2,WFAligner::Alignment,WFAligner::MemoryHigh);

    string pattern = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT";
    string text    = "TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
    aligner.alignEnd2End(pattern,text); // Align

    // Display CIGAR & score
    string cigar = aligner.getAlignmentCigar();
    cout << "CIGAR: " << cigar  << endl;
    cout << "Alignment score " << aligner.getAlignmentScore() << endl;
    return 0;
}
