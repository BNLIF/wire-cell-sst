#ifndef WIRECELLSST_DETECTORGDS
#define WIRECELLSST_DETECTORGDS

#include "WCPNav/DetectorGDS.h"
#include "WCPNav/GeomDataSource.h"
#include "WCPNav/WrappedGDS.h"

#include <istream>
#include <string>
#include <vector>

namespace WCPSst {

    /**
       WCPSst::DetectorGDS - read in a Channel Wire Geometry database.
       
     */

    class DetectorGDS : public WCP::DetectorGDS {//WCP::DetectorGDS {
    public:
	/// Read from an input stream containing content from ChannelWireGeometry.txt
        DetectorGDS(std::vector<std::string> geometry);
	virtual ~DetectorGDS();

    };

} 

#endif
