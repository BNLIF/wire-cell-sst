#include "WCPSst/GeomDataSource.h"
#include <iostream>
using namespace WCP;
using namespace std;

int main(int argc, char* argv[])
{
    WCPSst::GeomDataSource gds(argv[1]);

    std::vector<double> ex = gds.extent();
    cerr << "Extent: "
	 << " x:" << ex[0]/units::mm << " mm"
	 << " y:" << ex[1]/units::m << " m"
	 << " z:" << ex[2]/units::m << " m"
	 << endl;

    return 0;
}
