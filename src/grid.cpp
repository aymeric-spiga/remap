#include "grid.h"

Coord readPole(std::istream &is)
{
	double lon, lat;
	is >> lon >> lat;
	return (lon<=380 && lat<=380) ? xyz(lon, lat) : ORIGIN;
}

CGrid srcGrid;
CGrid tgtGrid;
