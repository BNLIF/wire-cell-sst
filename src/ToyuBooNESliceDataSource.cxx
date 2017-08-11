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
    , gds_flag(0)
    , gds(0)
    , nwire_u(0)
    , nwire_v(0)
    , nwire_w(0)
    , threshold_u(0)
    , threshold_v(0)
    , threshold_w(0)
    , threshold_ug(0)
    , threshold_vg(0)
    , threshold_wg(0)
{
    update_slices_bounds();
}

void WireCellSst::ToyuBooNESliceDataSource::set_flag(int flag1){
  flag = flag1;
}

WireCellSst::ToyuBooNESliceDataSource::ToyuBooNESliceDataSource(FrameDataSource& fds, FrameDataSource& fds1, float th_u, float th_v, float th_w, float th_ug, float th_vg, float th_wg, int nwire_u, int nwire_v, int nwire_w, std::vector<float>* uplane_rms, std::vector<float>* vplane_rms, std::vector<float>* wplane_rms)
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
  , uplane_rms(uplane_rms)
  , vplane_rms(vplane_rms)
  , wplane_rms(wplane_rms)
  , gds_flag(0)
  , gds(0)
{
    update_slices_bounds();
}


WireCellSst::ToyuBooNESliceDataSource::ToyuBooNESliceDataSource(const DetectorGDS& gds, FrameDataSource& fds, FrameDataSource& fds1, float th_u, float th_v, float th_w, float th_ug, float th_vg, float th_wg, int nwire_u, int nwire_v, int nwire_w, std::vector<float>* uplane_rms, std::vector<float>* vplane_rms, std::vector<float>* wplane_rms)
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
  , uplane_rms(uplane_rms)
  , vplane_rms(vplane_rms)
  , wplane_rms(wplane_rms)
  , gds_flag(1)
  , gds(&gds)
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
	saved_signal[trace.chid] = q;
	
	if (q>threshold){
	  slice_group.push_back(Channel::Charge(trace.chid, q));
	  fired_channels.insert(trace.chid);
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
    }else if (flag==1){
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
	int q_next=0, q_prev=0;
	int q1 = trace1.charge[slice_tbin];
	// int q1_next, q1_prev;
	saved_signal[trace.chid] = q1;

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

	float threshold_g=0;

	if (uplane_rms ==0 && vplane_rms ==0 && wplane_rms == 0 ){

	  if (gds_flag == 0){
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
	  }else{
	    WirePlaneType_t plane = gds->channel_plane_conv(trace.chid);
	    if (plane == WirePlaneType_t(0)){
	      threshold =threshold_u;
	      threshold_g = threshold_ug;
	    }else if (plane == WirePlaneType_t(1)){
	      threshold = threshold_v;
	      threshold_g = threshold_vg;
	    }else if (plane == WirePlaneType_t(2)){
	      threshold = threshold_w;
	      threshold_g = threshold_wg;
	    }
	  }

	}else{
	  // Note: we did not consider the rebin here, so the threshold is effectively low ... 
	  if (gds_flag == 0){
	    if (trace.chid < nwire_u){
	      threshold = 3.6 * (*uplane_rms).at(trace.chid); // 3.6 sigma?
	      threshold_g = threshold_ug;
	      if (threshold == 0 ) threshold = threshold_u;
	    }else if (trace.chid < nwire_u + nwire_v){
	      threshold = 3.6 * (*vplane_rms).at(trace.chid - nwire_u);
	      threshold_g = threshold_vg;
	      if (threshold == 0 ) threshold = threshold_v;
	    }else if (trace.chid < nwire_u + nwire_v + nwire_w){
	      threshold = 3.6 * (*wplane_rms).at(trace.chid - nwire_u - nwire_v);
	      threshold_g = threshold_wg;
	      if (threshold == 0 ) threshold = threshold_w;
	    }
	  }else{
	    WirePlaneType_t plane = gds->channel_plane_conv(trace.chid);
	    if (plane == WirePlaneType_t(0)){
	      threshold =threshold_u;
	      threshold_g = threshold_ug;
	    }else if (plane == WirePlaneType_t(1)){
	      threshold = threshold_v;
	      threshold_g = threshold_vg;
	    }else if (plane == WirePlaneType_t(2)){
	      threshold = threshold_w;
	      threshold_g = threshold_wg;
	    }
	  }

	}

	//std::cout << q << " " << q1 << " " << threshold << std::endl;

	if (q>threshold){
	  slice_group.push_back(Channel::Charge(trace.chid, q1));
	  fired_channels.insert(trace.chid);
	  //	}else{
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
	  //	}
	}else if (q_next > threshold || q_prev > threshold){
	  if ((q1 > q_next/3.&&q_next>threshold) || (q1 > q_prev/3.&&q_prev>threshold)) {// scale by a factor of 3 ... 
	    slice_group.push_back(Channel::Charge(trace.chid, q1));
	    fired_channels.insert(trace.chid);
	  }
	}

	// else if (q_next > threshold || q_prev > threshold){
	//   if (q1 > threshold_g){
	//     slice_group.push_back(Channel::Charge(trace.chid, q1));
	//   }
	// }

	//std::cout << trace.chid << " " << q/threshold << " " << q1 << std::endl;
	
      }
    }else if (flag==2){
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
	
	WirePlaneType_t plane = gds->channel_plane_conv(trace.chid);
	if (plane == WirePlaneType_t(0)){
	  threshold =threshold_u;
	}else if (plane == WirePlaneType_t(1)){
	  threshold = threshold_v;
	}else if (plane == WirePlaneType_t(2)){
	  threshold = threshold_w;
	}

	// Save association of a channel ID and its charge.
	int q = trace.charge[slice_tbin];
	
	if (q>threshold){
	  //	  std::cout << q << " " << threshold << std::endl;
	  slice_group.push_back(Channel::Charge(trace.chid, q));
	  fired_channels.insert(trace.chid);
	}
      }
    }
    
    // if (flag!=2){
    //   // try to save stuff for the coherent subtraction ... 
    //   // compare the 48 channels
    //   // require >30 channels to fire?
    //   // require an up to 5 empty channels?
    //   // require both end has fired channels
    //   int step_limit = 5;
    //   int counter[172]={0};
    //   // for (int i=0;i!=172;i++){
    //   // 	counter[i] = 0;
    //   // }
    //   for (auto it = fired_channels.begin();it!=fired_channels.end();it++){
    // 	int count = int(*it/48);
    // 	counter[count]++;
    //   }
    //   for (int i=0;i!=172;i++){
    // 	if (counter[i]>30){
    // 	  //std::cout << "Xin " << slice_tbin << " " << i << " " << counter[i] << std::endl;
    // 	  for (int j=i*48;j!=(i+1)*48;j++){
    // 	    // this channel need to be not fired
    // 	    auto it = fired_channels.find(j);
    // 	    if (it == fired_channels.end()){
    // 	      int prev_channel,next_channel;
    // 	      // find the previous fired channel
    // 	      for (int k=j-1;k!=j-1-step_limit;k--){
    // 		prev_channel = k;
    // 		auto it1 = fired_channels.find(k);
    // 		if (it1 != fired_channels.end()) break;
    // 	      }
    // 	      // find the next fired channel
    // 	      for (int k=j+1;k!=j+1+step_limit;k++){
    // 		next_channel = k;
    // 		auto it1 = fired_channels.find(k);
    // 		if (it1 != fired_channels.end()) break;
    // 	      }
	      
    // 	      if (fabs(next_channel - prev_channel) <= step_limit+1){
    // 		slice_group.push_back(Channel::Charge(j, saved_signal[j]));
    // 	      }
    // 	    }
    // 	  }
	  
    // 	}
    //   }
    //   //
    // }


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

