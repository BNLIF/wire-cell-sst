#ifndef WIRECELLSST_UTIL
#define WIRECELLSST_UTIL

#include "WCPSst/FrameDataSource.h"

#include "TFile.h"

namespace WCPSst {

    /** Make and return a frame data source attached to the given filename.
     *
     * This will attempt to check for file schema versions.
     */
    WCPSst::FrameDataSource* make_fds(TFile& tfile, const char* tpath = "/Event/Sim");
    WCPSst::FrameDataSource* make_fds(const char* filename, const char* tpath = "/Event/Sim");


}

#endif

