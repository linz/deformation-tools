#ifndef GRIDUTIL_HPP
#define GRIDUTIL_HPP

#include "grid.hpp"

void mark_file( grid &g, const char *markfile );
void outline_affected_cells( grid &g, const char *wktfile );
// write_linz_grid columns is blank or list of columns to export separated by "+"
void write_linz_grid_file( grid &g, std::string crdsys, std::string header1, std::string header2, std::string header3, std::string vres, std::string columns, const char *gridfile );
void print_grid_stats( grid &g );
void write_grid_extents_wkt( grid &g, const char *wkt_file );
#endif
