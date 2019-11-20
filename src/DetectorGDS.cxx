#include "WCPSst/DetectorGDS.h"
#include "WCPData/Units.h"
#include "WCPData/GeomWire.h"

#include <string>
#include <sstream>
#include <cassert>
#include <fstream>
#include <vector>

using namespace units;


WCPSst::DetectorGDS::DetectorGDS(std::vector<std::string> geometry)
    : WCP::DetectorGDS(geometry)
{
  /*
    WCP::DetectorGDS det_gds(geometry);
    if (!geometry.empty()) {
        for (short cryo = 0; cryo < det_gds.ncryos(); cryo++) {
	    for (short apa = 0; apa < det_gds.napa(cryo); apa++) {
	      WCP::WrappedGDS *apa_gds = det_gds.get_apaGDS(cryo, apa);
	      std::cout<<"got wrappedGDS"<<std::endl;
	      this->load_apa((const WCP::WrappedGDS &)apa_gds);
	      std::cout<<"loaded wires"<<std::endl;
	    }
	}
    }
  */
}
/*
void WCPSst::DetectorGDS::load_apa(const WCP::WrappedGDS &apa_gds)
{

    for (int iplane = 0; iplane < 3; iplane++) {
        WCP::WirePlaneType_t plane = (WCP::WirePlaneType_t)iplane;
	std::cout<<"before getting wires in plane"<<std::endl;
	const WCP::GeomWireSelection &wip = apa_gds.wires_in_plane(plane);
	std::cout<<"got GeomWireSelection "<<wip.size()<<std::endl;
	for (auto wit = wip.begin(); wit != wip.end(); ++wit) {
	    this->add_wire((const WCP::GeomWire &)wit);
	}
    }
}
*/	
WCPSst::DetectorGDS::~DetectorGDS()
{
}

