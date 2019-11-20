#ifndef WIRECELLSST_UBOONESLICEDATASOURCE_H
#define WIRECELLSST_UBOONESLICEDATASOURCE_H

#include "WCPNav/FrameDataSource.h"
#include "WCPNav/GeomDataSource.h"
#include "WCPNav/DetectorGDS.h"
#include "WCPData/Slice.h"
#include "WCPData/Frame.h"
#include "WCPData/GeomWire.h"

namespace WCPSst {

    /**
       SliceDataSource - deliver slices of frames from a FrameDataSource
     */
    class uBooNESliceDataSource {
    public:

      uBooNESliceDataSource(WCP::FrameDataSource& fds, WCP::FrameDataSource& fds1, WCP::FrameDataSource& fds1_error, float th_u, float th_v, float th_w, int nwire_u, int nwire_v, int nwire_w, std::vector<float>* uplane_rms = 0, std::vector<float>* vplane_rms = 0, std::vector<float>* wplane_rms = 0);
      
      virtual ~uBooNESliceDataSource();
      
      /// Return the number of slices in the current frame.  
      virtual int size() const;
      /// Go to the given slice, return slice number or -1 on error
      virtual int jump(int slice_number); 
      
      /// Go to the next slice, return its number or -1 on error
      virtual int next();
      
      /// Get the current slice
      virtual WCP::Slice&  get();
      virtual const WCP::Slice&  get() const;

      virtual WCP::Slice& get_error();
      virtual const WCP::Slice&  get_error() const;
      
    private:
      
      WCP::FrameDataSource& _fds;
      WCP::FrameDataSource& _fds1;
      WCP::FrameDataSource& _fds2;
      
     
      int nwire_u, nwire_v, nwire_w;
      
      WCP::Slice _slice;	// cache the current slice
      WCP::Slice _slice_error;	// cache the current slice
      
      int _frame_index;	// last frame we loaded
      int _slice_index;	// current slice, for caching
      mutable int _slices_begin; // tbin index of earliest bin of all traces
      mutable int _slices_end; // tbin index of one past the latest bin of all traces
      float threshold;
      
      float threshold_u;
      float threshold_v;
      float threshold_w;
      
      std::vector<float>* uplane_rms;
      std::vector<float>* vplane_rms;
      std::vector<float>* wplane_rms;
      
      
      virtual void update_slices_bounds() const;
      virtual void clear();
      
    };
    
}

#endif
