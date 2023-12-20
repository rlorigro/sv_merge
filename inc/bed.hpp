#pragma once

#include "Region.hpp"

#include <vector>

using std::vector;


namespace sv_merge {

void for_region_in_bed_file(path bed_path, const function<void(const Region& r)>& f);

void load_windows_from_bed(path windows_bed, vector<Region>& regions);

}
