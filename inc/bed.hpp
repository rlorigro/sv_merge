#pragma once

#include "Region.hpp"


namespace sv_merge {

void for_region_in_bed_file(path bed_path, const function<void(const Region& r)>& f);

}