#ifndef WIRECELLSST_DETECTORGDS
#define WIRECELLSST_DETECTORGDS

#include "WireCellNav/DetectorGDS.h"
#include "WireCellNav/GeomDataSource.h"
#include "WireCellNav/WrappedGDS.h"

#include <istream>
#include <string>
#include <vector>

namespace WireCellSst {

    /**
       WireCellSst::DetectorGDS - read in a Channel Wire Geometry database.
       
     */

    class DetectorGDS : public WireCell::DetectorGDS {//WireCell::DetectorGDS {
    public:
	/// Read from an input stream containing content from ChannelWireGeometry.txt
        DetectorGDS(std::vector<std::string> geometry);
	virtual ~DetectorGDS();

    };

} 

#endif
