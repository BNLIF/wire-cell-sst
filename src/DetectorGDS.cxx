#include "WireCellSst/DetectorGDS.h"
#include "WireCellData/Units.h"
#include "WireCellData/GeomWire.h"

#include <string>
#include <sstream>
#include <cassert>
#include <fstream>
#include <vector>

using namespace units;


WireCellSst::DetectorGDS::DetectorGDS(std::vector<std::string> geometry)
    : WireCell::DetectorGDS(geometry)
{
  /*
    WireCell::DetectorGDS det_gds(geometry);
    if (!geometry.empty()) {
        for (short cryo = 0; cryo < det_gds.ncryos(); cryo++) {
	    for (short apa = 0; apa < det_gds.napa(cryo); apa++) {
	      WireCell::WrappedGDS *apa_gds = det_gds.get_apaGDS(cryo, apa);
	      std::cout<<"got wrappedGDS"<<std::endl;
	      this->load_apa((const WireCell::WrappedGDS &)apa_gds);
	      std::cout<<"loaded wires"<<std::endl;
	    }
	}
    }
  */
}
/*
void WireCellSst::DetectorGDS::load_apa(const WireCell::WrappedGDS &apa_gds)
{

    for (int iplane = 0; iplane < 3; iplane++) {
        WireCell::WirePlaneType_t plane = (WireCell::WirePlaneType_t)iplane;
	std::cout<<"before getting wires in plane"<<std::endl;
	const WireCell::GeomWireSelection &wip = apa_gds.wires_in_plane(plane);
	std::cout<<"got GeomWireSelection "<<wip.size()<<std::endl;
	for (auto wit = wip.begin(); wit != wip.end(); ++wit) {
	    this->add_wire((const WireCell::GeomWire &)wit);
	}
    }
}
*/	
WireCellSst::DetectorGDS::~DetectorGDS()
{
}

