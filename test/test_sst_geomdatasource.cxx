#include "WCPSst/GeomDataSource.h"

#include <iostream>
#include <sstream>
#include <string>
#include <cmath>

using namespace WCP;
using namespace std;

int main()
{
    // R"()" is C++-0x
    // The tests below bake-in various things about this test data so don't go changing it willy-nilly.
    std::string data = R"(
# This file is used to translate between channel ID, wire ID, and geometry.
#
# channel: channel ID
# plane: plane number (0 == U, 1 == V, 2 == Z)
# wire: wire ID on that plane
# sx, sy, sz: starting (x,y,z) position of the wire in cm
# ex, ey, ez: ending (x,y,z) position of the wire in cm
#
# ch pl wire       sx     sy    sz          ex     ey   ez
0    0 0          0.0  117.021 0.206489    0.0  117.36 0.793511
1    0 1          0.0  116.675 0.206489    0.0  117.36 1.39351
2    0 2          0.0  116.328 0.206489    0.0  117.36 1.99351
3    0 3          0.0  115.982 0.206489    0.0  117.36 2.59351
4    0 4          0.0  115.636 0.206489    0.0  117.36 3.19351
5    0 5          0.0  115.289 0.206489    0.0  117.36 3.79351
6    0 6          0.0  114.943 0.206489    0.0  117.36 4.39351
7    0 7          0.0  114.596 0.206489    0.0  117.36 4.99351
8    0 8          0.0  114.25  0.206489    0.0  117.36 5.59351
9    0 9          0.0  113.903 0.206489    0.0  117.36 6.19351 
2399 1 0         -0.3 -115.081 0.00648934 -0.3 -115.42 0.593511
2400 1 1         -0.3 -114.735 0.00648934 -0.3 -115.42 1.19351
2401 1 2         -0.3 -114.388 0.00648934 -0.3 -115.42 1.79351
2402 1 3         -0.3 -114.042 0.00648934 -0.3 -115.42 2.39351
2403 1 4         -0.3 -113.696 0.00648934 -0.3 -115.42 2.99351
2404 1 5         -0.3 -113.349 0.00648934 -0.3 -115.42 3.59351
2405 1 6         -0.3 -113.003 0.00648934 -0.3 -115.42 4.19351
2406 1 7         -0.3 -112.656 0.00648934 -0.3 -115.42 4.79351
2407 1 8         -0.3 -112.31  0.00648934 -0.3 -115.42 5.39351
2408 1 9         -0.3 -111.963 0.00648934 -0.3 -115.42 5.99351 
4798 2 0         -0.6 -115.53  0.3        -0.6  117.47 0.3
4799 2 1         -0.6 -115.53  0.6        -0.6  117.47 0.6
4800 2 2         -0.6 -115.53  0.9        -0.6  117.47 0.9
4801 2 3         -0.6 -115.53  1.2        -0.6  117.47 1.2
4802 2 4         -0.6 -115.53  1.5        -0.6  117.47 1.5
4803 2 5         -0.6 -115.53  1.8        -0.6  117.47 1.8
4804 2 6         -0.6 -115.53  2.1        -0.6  117.47 2.1
4805 2 7         -0.6 -115.53  2.4        -0.6  117.47 2.4
4806 2 8         -0.6 -115.53  2.6        -0.6  117.47 2.8
4807 2 9         -0.6 -115.53  3          -0.6  117.47 3)";

    istringstream iss(data);
    WCPSst::GeomDataSource gds;
    gds.load(iss);

    int errors = 0;
    
    for (int iplane = 0; iplane < 3; ++iplane) {
	WirePlaneType_t plane = static_cast<WirePlaneType_t>(iplane);
	GeomWireSelection ws = gds.wires_in_plane(plane);
	if (ws.size() != 10) {
	    cerr << "FAIL: got wrong number of wires." << endl;
	    ++errors;
	}
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
		 << " plane=" << wire.plane() << " index=" << wire.index() << endl;
	}

	double pitch = gds.pitch(plane);
	double angle = gds.angle(plane);
	cout << "pitch=" << pitch/units::mm << " mm"
	     <<" angle=" << angle/units::degree << " degree"
	     << endl;
	if (std::isnan(pitch)) {
	    cerr << "FAIL: got NaN for pitch for plane " << plane << endl;
	    ++errors;
	}
	if (std::abs(pitch-3*units::mm) > 0.01*units::mm) {
	    cerr << "FAIL: got wrong pitch for plane " << plane << ": " << pitch << endl;
	    ++errors;
	}
	double wantang[3] = {60.0*units::degree, -60.0*units::degree, 0.0*units::degree};
	double want = wantang[iplane];
	double epsdegree = 0.01*units::degree;
	if (std::abs(angle - want) > epsdegree) {
	    cerr << "FAIL: got wrong angle for plane " << plane << ": " << angle/units::degree << " degree"
		 << " wanted: " << want/units::degree
		 << endl;
	    ++errors;
	}

	for (int mmind=0; mmind<3; ++mmind) {
	    std::pair<double,double> mm = gds.minmax(mmind, plane);
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


    }


    for (int iplane=-1; iplane<3; ++iplane) {
        std::vector<double> ex = gds.extent((WCP::WirePlaneType_t)iplane);
        cerr << "Extent for plane "<<iplane<<": "
             << " x:" << ex[0]/units::mm << " mm,"
             << " y:" << ex[1]/units::mm << " mm,"
             << " z:" << ex[2]/units::mm << " mm"
             << endl;
    }

    std::vector<double> ex = gds.extent(WCP::kUnknownWirePlaneType);
    cerr << "Extent: "
	 << " x:" << ex[0]/units::mm << " mm,"
	 << " y:" << ex[1]/units::mm << " mm,"
	 << " z:" << ex[2]/units::mm << " mm"
	 << endl;
    // note we are using only a portion of wires, which explains these weird numbers
    double wantext[] = {6.0*units::mm, 2.33*units::m, 6.19*units::cm};
    for (int ind=0; ind<3; ++ind) {
        std::pair<double,double> mm = gds.minmax(ind);

	double normed = std::abs(ex[ind] + wantext[ind]);
        if (normed) {
            normed = std::abs(ex[ind] - wantext[ind]) / normed;
        }
	cerr << "Extent for all planes: dim#" << ind << " got-want: " << ex[ind] << " - " << wantext[ind]
             << " normed: " << normed
             << " totext=[" << mm.first/units::mm << "mm , " << mm.second/units::mm << "mm]" 
             << endl;
	if (normed > 0.001) {
	    cerr << "FAIL: got unexpected extent for dimension " << ind
                 << ", plane all: " << normed << endl;
	    ++errors;
	}

    }



    if (errors) {
	cerr << "FAILED: got " << errors << " errors" << endl;
	exit(1);
    }

    return 0;
}
