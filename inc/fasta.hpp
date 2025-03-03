#pragma once

#include <functional>
#include <filesystem>
using std::filesystem::path;
using std::function;

#include "Sequence.hpp"

namespace sv_merge {

void for_sequence_in_fasta_file(path fasta_path, const function<void(const Sequence& s)>& f);

} // sv_merge
