#include "WireCellSst/ToyuBooNESliceDataSource.h"

using namespace WireCell;

WireCellSst::ToyuBooNESliceDataSource::ToyuBooNESliceDataSource(FrameDataSource& fds, float th)
    : _fds(fds)
    , _fds1(fds)
    , _frame_index(-1)
    , _slice_index(-1)
    , _slices_begin(-1)
    , _slices_end(-1)
    , threshold(th)
    , flag(0)
{
    update_slices_bounds();
}

WireCellSst::ToyuBooNESliceDataSource::ToyuBooNESliceDataSource(FrameDataSource& fds, FrameDataSource& fds1, float th_u, float th_v, float th_w, float th_ug, float th_vg, float th_wg, int nwire_u, int nwire_v, int nwire_w)
    : _fds(fds)
    , _fds1(fds1)
    , _frame_index(-1)
    , _slice_index(-1)
    , _slices_begin(-1)
    , _slices_end(-1)
    , threshold_u(th_u)
    , threshold_v(th_v)
    , threshold_w(th_w)
    , threshold_ug(th_ug)
    , threshold_vg(th_vg)
    , threshold_wg(th_wg)
    , flag(1)
    , nwire_u(nwire_u)
    , nwire_v(nwire_v)
    , nwire_w(nwire_w)
{
    update_slices_bounds();
}


WireCellSst::ToyuBooNESliceDataSource::~ToyuBooNESliceDataSource()
{
}

void WireCellSst::ToyuBooNESliceDataSource::clear()
{
    _slice_index = _slices_begin = _slices_end = -1;
    _slice.clear();
}

// Update the bounds on the slice indices based on the current frame.
void WireCellSst::ToyuBooNESliceDataSource::update_slices_bounds() const
{
    const Frame& frame = _fds.get();

    if (frame.index < 0) {	// no frames
	return;
    }

    if (frame.index == _frame_index) { // no change, no-op
	return;
    }

    // frame change

    size_t ntraces = frame.traces.size();
    for (size_t ind=0; ind<ntraces; ++ind) {
	const Trace& trace = frame.traces[ind];

	if (!ind) {		// first time
	    _slices_begin = trace.tbin;
	    _slices_end   = trace.charge.size() + _slices_begin;
	    continue;
	}
	int tbin = trace.tbin;
	int nbins = trace.charge.size();
	_slices_begin = std::min(_slices_begin, tbin);
	_slices_end   = std::max(_slices_end,   tbin + nbins);
    }
    return;
}

int WireCellSst::ToyuBooNESliceDataSource::size() const
{
    this->update_slices_bounds();
    if (_slices_begin < 0 || _slices_end < 0) {
	return 0;
    }
    return _slices_end - _slices_begin;
}

int WireCellSst::ToyuBooNESliceDataSource::jump(int index)
{
    if (index < 0) {		// underflow
	this->clear();
	return index;
    }

    this->update_slices_bounds();

    if (index >= _slices_end) {	// overflow
	this->clear();
	return -1;
    }
    if (index == _slice_index) {
	return index;		// already loaded, no-op
    }

    // new slice

    _slice.clear();
    _slice_index = index;	
    int slice_tbin = _slice_index + _slices_begin; 
    Channel::Group slice_group;
    
    //int sum = 0;

    if (flag==0){
      const Frame& frame = _fds.get();
      size_t ntraces = frame.traces.size();
      for (size_t ind=0; ind<ntraces; ++ind) {
	const Trace& trace = frame.traces[ind];
	int tbin = trace.tbin;
	int nbins = trace.charge.size();
	
	if (slice_tbin < tbin) {
	  continue;
	}
	if (slice_tbin >= tbin+nbins) {
	  continue;
	}
	
	// Save association of a channel ID and its charge.
	int q = trace.charge[slice_tbin];
	if (q>threshold){
	  slice_group.push_back(Channel::Charge(trace.chid, q));
	}
	// else{
	//   if (umap!=0&&vmap!=0&&wmap!=0){
	//     // hack for now for data
	//     int nwire_u = 2400;
	//     int nwire_v = 2400;
	//     int nwire_w = 3456;
	//     if (trace.chid < nwire_u){
	//       if (umap->find(trace.chid) == umap->end())
	//   	slice_group.push_back(Channel::Charge(trace.chid, q));
	//     }else if (trace.chid < nwire_u + nwire_v){
	//       if (vmap->find(trace.chid-nwire_u) == vmap->end())
	//   	slice_group.push_back(Channel::Charge(trace.chid, q));
	//     }else{
	//       if (wmap->find(trace.chid-nwire_u-nwire_v) == wmap->end())
	//   	slice_group.push_back(Channel::Charge(trace.chid, q));
	//     }
	//   }
	  // hack for now
	//	}

      }
    }else{
      // may need to update to take into account that the two frames may not be synced ... 
      const Frame& frame = _fds.get();
      const Frame& frame1 = _fds1.get();
      size_t ntraces = frame.traces.size();
      for (size_t ind=0; ind<ntraces; ++ind) {
	const Trace& trace = frame.traces[ind];
	const Trace& trace1 = frame1.traces[ind];
	
	int tbin = trace.tbin;
	int nbins = trace.charge.size();
	
	if (slice_tbin < tbin) {
	  continue;
	}
	if (slice_tbin >= tbin+nbins) {
	  continue;
	}
	
	// Save association of a channel ID and its charge.
	int q = trace.charge[slice_tbin];
	int q_next, q_prev;
	int q1 = trace1.charge[slice_tbin];
	// int q1_next, q1_prev;
	
	if (slice_tbin == tbin){
	  q_next  = trace.charge[slice_tbin +1];
	  q_prev  = trace.charge[slice_tbin +1];
	  // q1_next = trace1.charge[slice_tbin +1];
	  // q1_prev = trace1.charge[slice_tbin +1];
	}else if (slice_tbin == tbin + nbins -1){
	  q_next  = trace.charge[slice_tbin -1];
	  q_prev  = trace.charge[slice_tbin -1];
	  // q1_next = trace1.charge[slice_tbin -1];
	  // q1_prev = trace1.charge[slice_tbin -1];
	}else{
	  q_next  = trace.charge[slice_tbin +1];
	  q_prev  = trace.charge[slice_tbin -1];
	  // q1_next = trace1.charge[slice_tbin +1];
	  // q1_prev = trace1.charge[slice_tbin -1];
	}

	float threshold_g;
	if (trace.chid < nwire_u){
	  threshold =threshold_u;
	  threshold_g = threshold_ug;
	}else if (trace.chid < nwire_u + nwire_v){
	  threshold = threshold_v;
	  threshold_g = threshold_vg;
	}else if (trace.chid < nwire_u + nwire_v + nwire_w){
	  threshold = threshold_w;
	  threshold_g = threshold_wg;
	}

	if (q>threshold){
	  slice_group.push_back(Channel::Charge(trace.chid, q1));
	}else{
	  // if (umap!=0&&vmap!=0&&wmap!=0){
	  //   // hack for now for data
	  //   int nwire_u = 2400;
	  //   int nwire_v = 2400;
	  //   int nwire_w = 3456;
	  //   if (trace.chid < nwire_u){
	  //     if (umap->find(trace.chid) == umap->end())
	  // 	slice_group.push_back(Channel::Charge(trace.chid, 0));
	  //   }else if (trace.chid < nwire_u + nwire_v){
	  //     if (vmap->find(trace.chid-nwire_u) == vmap->end())
	  // 	slice_group.push_back(Channel::Charge(trace.chid, 0));
	  //   }else{
	  //     if (wmap->find(trace.chid-nwire_u-nwire_v) == wmap->end())
	  // 	slice_group.push_back(Channel::Charge(trace.chid, 0));
	  //   }
	  // }
	  //   //hack for now 
	}


// else if (q_next > threshold || q_prev > threshold){
	//   if (q1 > threshold_g){
	//     slice_group.push_back(Channel::Charge(trace.chid, q1));
	//   }
	// }

	//std::cout << trace.chid << " " << q/threshold << " " << q1 << std::endl;
	
      }
    }
    
    _slice.reset(slice_tbin, slice_group);
    //std::cout << sum << std::endl;
    return index;
}

int WireCellSst::ToyuBooNESliceDataSource::next()
{
    this->jump(_slice_index+1);
}

Slice& WireCellSst::ToyuBooNESliceDataSource::get()
{
    return _slice;
}

const Slice& WireCellSst::ToyuBooNESliceDataSource::get() const
{
    return _slice;
}

