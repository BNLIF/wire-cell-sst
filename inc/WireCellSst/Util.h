#ifndef WIRECELLSST_UTIL
#define WIRECELLSST_UTIL

#include "WireCellSst/FrameDataSource.h"

#include "TFile.h"

namespace WireCellSst {

    /** Make and return a frame data source attached to the given filename.
     *
     * This will attempt to check for file schema versions.
     */
    WireCellSst::FrameDataSource* make_fds(TFile& tfile, const char* tpath = "/Event/Sim");
    WireCellSst::FrameDataSource* make_fds(const char* filename, const char* tpath = "/Event/Sim");


}

#endif

