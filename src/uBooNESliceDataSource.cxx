#include "WireCellSst/uBooNESliceDataSource.h"

using namespace WireCell;



WireCellSst::uBooNESliceDataSource::uBooNESliceDataSource(FrameDataSource& fds, FrameDataSource& fds1,  FrameDataSource& fds2, float th_u, float th_v, float th_w, int nwire_u, int nwire_v, int nwire_w, std::vector<float>* uplane_rms, std::vector<float>* vplane_rms, std::vector<float>* wplane_rms)
    : _fds(fds)
    , _fds1(fds1)
    , _fds2(fds2)
    , _frame_index(-1)
    , _slice_index(-1)
    , _slices_begin(-1)
    , _slices_end(-1)
    , threshold_u(th_u)
    , threshold_v(th_v)
    , threshold_w(th_w)
    , nwire_u(nwire_u)
    , nwire_v(nwire_v)
    , nwire_w(nwire_w)
    , uplane_rms(uplane_rms)
    , vplane_rms(vplane_rms)
    , wplane_rms(wplane_rms)
{
    update_slices_bounds();
}


WireCellSst::uBooNESliceDataSource::~uBooNESliceDataSource()
{
}

void WireCellSst::uBooNESliceDataSource::clear()
{
    _slice_index = _slices_begin = _slices_end = -1;
    _slice.clear();
}

// Update the bounds on the slice indices based on the current frame.
void WireCellSst::uBooNESliceDataSource::update_slices_bounds() const
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

int WireCellSst::uBooNESliceDataSource::size() const
{
    this->update_slices_bounds();
    if (_slices_begin < 0 || _slices_end < 0) {
	return 0;
    }
    return _slices_end - _slices_begin;
}

int WireCellSst::uBooNESliceDataSource::jump(int index)
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


    int saved_signal[8256]={0};
    for (int i=0;i!=8256;i++){
      saved_signal[i] = 0;
    }
    std::set<int> fired_channels;
    // new slice

    _slice.clear();
    _slice_index = index;	
    int slice_tbin = _slice_index + _slices_begin; 
    Channel::Group slice_group;
    Channel::Group slice_group_error;
    
    //int sum = 0;

   
    // may need to update to take into account that the two frames may not be synced ... 
    const Frame& frame = _fds.get();
    const Frame& frame1 = _fds1.get();
    const Frame& frame2 = _fds2.get();
    
    size_t ntraces = frame.traces.size();
    for (size_t ind=0; ind<ntraces; ++ind) {
      const Trace& trace = frame.traces[ind];
      const Trace& trace1 = frame1.traces[ind];
      const Trace& trace2 = frame2.traces[ind];
      
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
      int q_next=0, q_prev=0;
      int q1 = trace1.charge[slice_tbin];
      int q1_error = trace2.charge[slice_tbin];
      // int q1_next, q1_prev;
      saved_signal[trace.chid] = q1;
      
      if (slice_tbin == tbin){
	q_next  = trace.charge[slice_tbin +1];
	q_prev  = trace.charge[slice_tbin +1];
      }else if (slice_tbin == tbin + nbins -1){
	q_next  = trace.charge[slice_tbin -1];
	q_prev  = trace.charge[slice_tbin -1];
      }else{
	q_next  = trace.charge[slice_tbin +1];
	q_prev  = trace.charge[slice_tbin -1];
      }
      
      if (trace.chid < nwire_u){
	threshold = 3.6 * (*uplane_rms).at(trace.chid); // 3.6 sigma?
	if (threshold == 0 ) threshold = threshold_u;
      }else if (trace.chid < nwire_u + nwire_v){
	threshold = 3.6 * (*vplane_rms).at(trace.chid - nwire_u);
	if (threshold == 0 ) threshold = threshold_v;
      }else if (trace.chid < nwire_u + nwire_v + nwire_w){
	threshold = 3.6 * (*wplane_rms).at(trace.chid - nwire_u - nwire_v);
	if (threshold == 0 ) threshold = threshold_w;
      }

      
      if (q>threshold){
	slice_group.push_back(Channel::Charge(trace.chid, q1));
	slice_group_error.push_back(Channel::Charge(trace.chid, q1_error));
	fired_channels.insert(trace.chid);
      }else if (q_next > threshold || q_prev > threshold){
	if ((q1 > q_next/3.&&q_next>threshold) || (q1 > q_prev/3.&&q_prev>threshold)) {// scale by a factor of 3 ... 
	  slice_group.push_back(Channel::Charge(trace.chid, q1));
	  slice_group_error.push_back(Channel::Charge(trace.chid, q1_error));
	  fired_channels.insert(trace.chid);
	}
      }
      
    }

    _slice.reset(slice_tbin, slice_group);
    _slice_error.reset(slice_tbin, slice_group_error);
    
    return index;
}

int WireCellSst::uBooNESliceDataSource::next()
{
    this->jump(_slice_index+1);
}

Slice& WireCellSst::uBooNESliceDataSource::get()
{
    return _slice;
}

const Slice& WireCellSst::uBooNESliceDataSource::get() const
{
    return _slice;
}


Slice& WireCellSst::uBooNESliceDataSource::get_error()
{
    return _slice_error;
}

const Slice& WireCellSst::uBooNESliceDataSource::get_error() const
{
    return _slice_error;
}

