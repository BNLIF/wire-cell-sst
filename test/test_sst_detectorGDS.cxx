#include "WireCellSst/DetectorGDS.h"
#include "WireCellNav/WrappedGDS.h"

#include <iostream>
#include <vector>
#include <string>

using namespace WireCell;
using namespace std;

int main()
{

    vector<string> geometry;
    geometry.push_back("/home/xiaoyueli/BNLIF/wire-cell/geom_35t_v5.txt");

    WireCellSst::DetectorGDS det_gds(geometry);
    
    int errors = 0;

    for (short cryo = 0; cryo < det_gds.ncryos(); cryo++) {
      for (short apa = 1; apa < 2/*det_gds.napa(cryo)*/; apa++) {
	  WrappedGDS *gds = det_gds.get_apaGDS(cryo, apa);
	  for  (int iplane = 0; iplane < 3; ++iplane) {
	      WirePlaneType_t plane = static_cast<WirePlaneType_t>(iplane);
	      GeomWireSelection ws = gds->wires_in_plane(plane);
	      sort_by_planeindex(ws);
	      for (size_t ind=0; ind<ws.size(); ++ind) {
		  const GeomWire& wire = *ws[ind];
		  if (wire.plane() != plane) {
		      cerr << "FAIL: got wrong plane" << endl;
		      ++errors;
		  }
		  if (wire.ident() < 0) {
		    cerr << "FAIL: bogus wire ID: " << wire 
			 << " plane=" << wire.plane() << " index=" << wire.index() << endl;
		    ++errors;
		  }
		  if (wire.index() < 0) {
		    cerr << "FAIL: bogus wire index: " << wire 
			 << " plane=" << wire.plane() << " index=" << wire.index() << endl;
		    ++errors;
		  }
		  
		  cout << "Wire: " << wire 
		    //<< " plane=" << wire.plane() << " index=" << wire.index() << endl;
		       << " x1 = " << (wire.point1()).x << ", x2 = " << (wire.point2()).x << endl;
	      }
	      
	      double pitch = gds->pitch(plane);
	      double angle = gds->angle(plane);
	      cout << "pitch=" << pitch/units::mm << " mm"
		   <<" angle=" << angle/units::degree << " degree"
		   << endl;
	      if (std::isnan(pitch)) {
		cerr << "FAIL: got NaN for pitch for plane " << plane << endl;
		++errors;
	      }
	      if (std::abs(pitch-4.88049*units::mm) > 0.01*units::mm) {
		cerr << "FAIL: got wrong pitch for plane " << plane << ": " << pitch << endl;
		++errors;
	      }
	      double wantang[3] = {-45.7*units::degree, 44.3*units::degree, 0.0*units::degree};
	      double want = wantang[iplane];
	      double epsdegree = 0.1*units::degree;
	      if (std::abs(angle - want) > epsdegree) {
		cerr << "FAIL: got wrong angle for plane " << plane << ": " << angle/units::degree << " degree"
		     << " wanted: " << want/units::degree
		     << endl;
		++errors;
	      }
	      
	      for (int mmind=0; mmind<3; ++mmind) {
		std::pair<double,double> mm = gds->minmax(mmind, plane);
		cerr << "plane=" << plane << "[" <<  mmind << "]"
		     << " min=" << mm.first/units::mm << " mm,"
		     << " max=" << mm.second/units::mm << " mm"
		     << endl;
		if (mmind && mm.first == mm.second && mm.first == 0.0) {
		  cerr << "FAIL: no mimmax!" << endl;
		  ++errors;
		}
		if (mm.first > mm.second) {
		  cerr << "FAIL: min bigger than max!" << endl;
		  ++errors;
		}
	      }	      	      

	      std::vector<double> ex = gds->extent((WireCell::WirePlaneType_t)iplane);
	      cerr << "Extent for plane "<<iplane<<": "
		   << " x:" << ex[0]/units::mm << " mm,"
		   << " y:" << ex[1]/units::mm << " mm,"
		   << " z:" << ex[2]/units::mm << " mm"
		   << endl;
	  }

	  std::vector<double> ex = gds->extent(WireCell::kUnknownWirePlaneType);
	  cerr << "Extent: "
	       << " x:" << ex[0]/units::mm << " mm,"
	       << " y:" << ex[1]/units::mm << " mm,"
	       << " z:" << ex[2]/units::mm << " mm"
	       << endl;

	}
    }

    if (errors) {
	cerr << "FAILED: got " << errors << " errors" << endl;
	exit(1);
    }

    return 0;
    
}
