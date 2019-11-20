#include "WCPSst/Util.h"

#include "TFile.h"
#include "TTree.h"

#include <string>


WCPSst::FrameDataSource* WCPSst::make_fds(const char* filename, const char* tpath)
{
    TFile* tfile = TFile::Open(filename, "READONLY");
    if (!tfile) {
	return 0;
    }
    return make_fds(*tfile, tpath);
}

WCPSst::FrameDataSource* WCPSst::make_fds(TFile& tfile, const char* tpath)
{
    TTree* sst = dynamic_cast<TTree*>(tfile.Get(tpath));
    if (!sst) {
	return 0;
    }

    //sst->SetDirectory(0);

    WCPSst::FrameDataSource* fds = new WCPSst::FrameDataSource(*sst);

    // now try to handle file schema version

    TObject* obj = tfile.Get("version");
    if (!obj) {
	return fds;
    }

    TNamed* named = dynamic_cast<TNamed*>(obj);
    if (!named) {
	return fds;
    }

    if (named->GetName() != std::string("version")) {
	return fds;		// probably superfluous
    }
    
    std::string version = named->GetTitle();
    if (version < "2.0") {
	return fds;
    }

    fds->set_sim_tree(*sst);	// SST uses same tree for both types of data
    return fds;
}


