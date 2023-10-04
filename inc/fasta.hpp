#pragma once

#include <functional>

using std::function;

#include "Filesystem.hpp"
#include "Sequence.hpp"

using ghc::filesystem::path;

namespace sv_merge {

void for_sequence_in_fasta_file(path fasta_path, const function<void(const Sequence& s)>& f);

} // sv_merge
