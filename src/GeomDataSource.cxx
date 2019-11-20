#include "WCPSst/GeomDataSource.h"
#include "WCPData/Units.h"

#include <string>
#include <sstream>
#include <cassert>
#include <fstream>

using namespace units;


WCPSst::GeomDataSource::GeomDataSource(const char* filename)
    : WCP::GeomDataSource()
{
    if (filename) {
	std::ifstream geotext(filename);
	this->load(geotext);
    }
}

void WCPSst::GeomDataSource::load(std::istream& geo)
{

    std::string line;
    while (std::getline(geo, line)) {
	if (! line.size()) {
	    continue;
	}
	if (line[0] == '#') {
	    continue;
	}

	std::istringstream iss(line);

	float sx, sy, sz, ex, ey, ez;
	int iplane, channel, index;

	iss >> channel >> iplane >> index
	    >> sx >> sy >> sz >> ex >> ey >> ez;
	assert (index >= 0);

	int ident = (iplane+1)*10000 + index;
	WCP::WirePlaneType_t plane = static_cast<WCP::WirePlaneType_t>(iplane);

	this->add_wire(WCP::GeomWire(ident, plane, index, channel,
					  WCP::Point(sx*cm,sy*cm,sz*cm),
					  WCP::Point(ex*cm,ey*cm,ez*cm)));
    }
}
	
WCPSst::GeomDataSource::~GeomDataSource()
{
}

