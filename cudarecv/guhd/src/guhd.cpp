
#include <uhd/types/tune_request.hpp>
#include <uhd/utils/thread_priority.hpp>
#include <uhd/utils/safe_call.hpp>
#include <uhd/usrp/multi_usrp.hpp>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <ctime>

#include <boost/format.hpp>
#include <boost/thread.hpp>
#include <boost/filesystem.hpp>

#include "guhd.h"
#include "buffer.h"

#define F_L1 (1575.42e6)
#define F_L2 (1227.60e6)
#define ANALOG_GAIN (50)
#define BW_L1 (2.0e6) // smallest BW supported by DBSRX2
#define BW_L2 (8.0e6)
#define SETUP_TIME (1.5)
#define NEW_FILE_CYCLE_MIN (10.0)

Guhd::Guhd(
    const boost::program_options::variables_map& vm,
/*        const std::vector<size_t>& ip,
        const std::vector<size_t>& ltwo,
        const std::string& clk,
        const std::string& sync,*/
        const bool& mimo
/*        const double& rate,
        const boost::filesystem::path& path*/):
	// Constructor initialize list to set up members which don't have default constructors
    _run(false),_overflow_message(true),
    _clock(vm["clk"] .as<std::string>()),
    _sync (vm["sync"].as<std::string>()), _mimo(mimo),
    _rate(vm["rate"] .as<double>()), _runN(0){

    const std::vector<size_t>& ip = vm["ip"].as< std::vector<size_t> >();

    _chan   .reserve(ip.size());
    _p_fbuff.reserve(ip.size());
    //_p_sock .reserve(ip.size());

    std::ostringstream oss;
    for(size_t i = 0 ; i < ip.size() ; i++){
        if (i) oss << ",";
        oss << "addr" << i << "=192.168.10." << ip[i];
        _chan.push_back(i);
    } std::clog << "[UHD] args = " << oss.str() << std::endl;
    std::clog << "[UHD] Setting usrp object & streamer" << std::endl;

    _usrp = uhd::usrp::multi_usrp::make(oss.str());
    UHD_ASSERT_THROW(_usrp->get_rx_num_channels() == _chan.size());

	// The params here define the sample sizes in data transfer and recording.
	//
	// CPU_FORMAT can be 8, 16, 24, or 32 bits and is meant to match 
	// the data size used for baseband processing.
	// OTW_FORMAT can be 8 or 16 bits and is the size of the samples being sent
	// over-the-wire; 8 bits is faster, 16 bits is higher resolution.
	//
	// Both data types are complex.
	// Per guhd.h, this will be 16 bits for both, unless SC8 is defined.
	// Probably best to keep CPU_FORMAT matching OTW_FORMAT for processing speed...
	// not sure that CPU_FORMAT of 32 bits provides any increase in resolution.
	// MEM_TYPE takes care of this automatically the way that guhd.h is setup, 
	// so don't change without good reason, and use MEM_TYPE in type definitions everywhere.
    uhd::stream_args_t stream_args(CPU_FORMAT,OTW_FORMAT);
    stream_args.channels = _chan;
    _rx_stream  = _usrp->get_rx_stream(stream_args);
    _spb        = _rx_stream->get_max_num_samps(); // This seems to be a fixed value based on USRP architecture 
    // (treat this call like magic since there is no documentation how it decides what is the max size)

    // Setting up L2 channels
    _b_ltwo.insert(_b_ltwo.begin(),_chan.size(),false);
    for (auto const& c:vm["ltwo"].as< std::vector<size_t> >()){
        UHD_ASSERT_THROW(c < _chan.size());
        _b_ltwo[c] = true;
        //std::clog << "[UHD] Channel #" << c << " will operate at F_L2 = 1227.60 MHz" << std::endl;
    }

    for (auto const& c:_chan){
        _p_fbuff.push_back(new FileBuff(vm,c,_b_ltwo[c],&_run));
        //_p_sock.push_back (new SockBuff(1111 + i,&_run));
    }
}

Guhd::~Guhd(){
    Stop();
    for(auto const& ptr:_p_fbuff) delete ptr;
    //for(auto const& ptr:_p_sock) delete ptr;
}

int
Guhd::Start(){
    if (_run) return -1;
    _run = true;

    //allocate buffers to receive with samples (one buffer per channel)
    std::clog << "[UHD] Allocating buffers: "<< _spb <<" samples/buffer" <<std::endl;

	// Allocate buffers. This is a vector of type MEM_TYPE of the size of the 
	// maximum possible number of samples per buffer (rx_streamer->get_max_num_samps)
	// (MEM_TYPE can be either 8 or 16 bits based on defines in guhd.h)
	//
	// This syntax is initializing the outer vector with as many elements as there are channels,
	// and each element being a vector of length samples-per-buffer
    _buffs = std::vector<std::vector<std::complex<MEM_TYPE> > > (
        _chan.size() , std::vector<std::complex<MEM_TYPE> >(_spb)
    );

    for(auto const& i:_chan){
        std::clog<<"\n[UHD] Configuring Channel #" << i << std::endl;
        _buff_ptrs.push_back(&_buffs[i].front());
        set_rf(i,_b_ltwo[i]);
    }

    // clock set
    set_clock();
    std::clog << "[UHD] " << "Using Device: " << _usrp->get_pp_string() << std::endl;

    uhd::stream_cmd_t stream_cmd( uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS);
    stream_cmd.num_samps = 0;
    stream_cmd.stream_now = false;
    stream_cmd.time_spec = uhd::time_spec_t(SETUP_TIME) + _usrp->get_time_now();
    _rx_stream->issue_stream_cmd(stream_cmd);

    auto t = std::time(nullptr);
    for (auto const &i:_chan){
        _p_fbuff[i]->/*init()->*/open(t);
        //_p_sock [i]->init()->open();
    }

    return 0;
}

int
Guhd::Update (){
    if (!_run) return -1;
    uhd::rx_metadata_t md;
    const size_t num_rx_samps = _rx_stream->recv(
        _buff_ptrs, _spb, md, 3.0 //timeout
    );

    //handle the error code
    if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_TIMEOUT) {
        std::cerr << "[UHD] " << boost::format("Timeout while streaming") << std::endl;
        return -1;
    }

    if (_overflow_message && md.error_code == uhd::rx_metadata_t::ERROR_CODE_OVERFLOW){
        _overflow_message = false;

        std::cerr << "[UHD] "
            << boost::format(
            "Got an overflow indication. Please consider the following:\n"
            "  Your write medium must sustain a rate of %fMB/s.\n"
            "  Dropped samples will not be written to the file.\n"
            "  This message will not appear again.\n"
        ) % (_usrp->get_rx_rate()*sizeof(std::complex<MEM_TYPE>) * _chan.size() /1e6);
    }

    if (md.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE){
        auto t = std::time(nullptr);
        auto tm = *std::gmtime(&t);

        const std::string error = str(boost::format("[UHD] Receiver error: %s") % md.strerror());
        std::cerr
            << std::put_time(&tm, "%Y%m%d %H:%M:%S ")
            << error << std::endl;
    }

	// _runN is the number of samples recorded in the file so far
	// So, if this reaches the minute-limit for the file, make a new file and reset
    if ((_runN += num_rx_samps) > _rate * 60 * NEW_FILE_CYCLE_MIN){
        auto t = std::time(nullptr);
        // If sample number limit has been exceeded, call the FileBuff member function
        // that will create a new file and update the pointer attribute for it
        for (auto const& i:_chan)
            _p_fbuff[i]->open(t);
        _runN = 0;
    }

	// No catastropic errors have occurred yet, so log the samples we do have
    for(auto const& i:_chan){
		// Call the FileBuff member function that will write this data to the current log file
        _p_fbuff[i]->append(_buffs[i],num_rx_samps);
    }

    return 0;
}

int
Guhd::Stop(){
    if (!_run) return -1;
    _run = false;
    for (auto const& i:_chan){
        _p_fbuff[i]/*->terminate()*/->close();
        //_p_sock[i] ->terminate()->close();
    }

    uhd::stream_cmd_t stream_cmd( uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS);
    _rx_stream->issue_stream_cmd(stream_cmd);
    std::clog << "[UHD] Stopping...";
    //for(auto& os : _ofs) os.close();
    //_ofs.clear();
    _buff_ptrs.clear();
    _buffs.clear();
    std::clog << " Stopped"<<std::endl;
    return 0;
}

void
Guhd::set_clock (){
    if (_mimo) {
        UHD_ASSERT_THROW(_usrp->get_num_mboards() == 2);
        std::clog << "[UHD] MIMO enabled: channel #1 set as slave." << std::endl;
        //make mboard 1 a slave over the MIMO Cable
        _usrp->set_clock_source("mimo", 1);
        _usrp->set_time_source ("mimo", 1);
    }

    size_t board_set = _mimo? 0 : uhd::usrp::multi_usrp::ALL_CHANS;
    _usrp->set_clock_source(_clock, board_set);

    if (_sync == "pps"){
        _usrp->set_time_source("external" , board_set );
        _usrp->set_time_unknown_pps(uhd::time_spec_t(0.0));
        boost::this_thread::sleep(boost::posix_time::seconds(1)); //wait for pps sync pulse
    } else //if (_sync == "internal"){
    //    _usrp->set_time_source_out(true,board_set);
        _usrp->set_time_now(uhd::time_spec_t(0.0));
    //} else UHD_ASSERT_THROW(false);


    // LO Lock Detect
    for (auto const& i:_chan){
        std::clog << "[UHD] Detecting clock lock: channel #" << i <<std::endl;
        check_locked_sensor(_usrp->get_rx_sensor_names(i), "lo_locked", boost::bind(&uhd::usrp::multi_usrp::get_rx_sensor, _usrp, _1, i));
        if (_mimo)
            check_locked_sensor(_usrp->get_mboard_sensor_names(i), "mimo_locked", boost::bind(&uhd::usrp::multi_usrp::get_mboard_sensor, _usrp, _1, i));
        if (_clock == "external")
            check_locked_sensor(_usrp->get_mboard_sensor_names(i), "ref_locked", boost::bind(&uhd::usrp::multi_usrp::get_mboard_sensor, _usrp, _1, i));
    }
}

void
Guhd::set_rf(const size_t& ch,const bool is_L2){
    uhd::tune_request_t tune_request(is_L2?F_L2:F_L1);
    _usrp->set_rx_freq(tune_request,ch);
    std::clog << "[UHD] " << boost::format("Actual RX Freq: %f MHz...") % (_usrp->get_rx_freq(ch)/1e6) << std::endl << std::endl;

    // SETTING RATE
    std::clog << "[UHD] " << boost::format("Setting RX Rate: %f Msps...") % (_rate/1e6) << std::endl;
    _usrp->set_rx_rate(_rate,ch);
    std::clog << "[UHD] " << boost::format("Actual RX Rate: %f Msps...") % (_usrp->get_rx_rate(ch)/1e6) << std::endl << std::endl;

    //set the rf gain
    std::clog << "[UHD] " << boost::format("Setting RX Gain: %f dB...") % ANALOG_GAIN << std::endl;
    _usrp->set_rx_gain(ANALOG_GAIN,ch);
    std::clog << "[UHD] " << boost::format("Actual RX Gain: %f dB...") % _usrp->get_rx_gain(ch) << std::endl << std::endl;

    //set the IF filter bandwidth
    const double bw = is_L2?BW_L2:BW_L1;
    std::clog << "[UHD] " << boost::format("Setting RX Bandwidth: %f MHz...") % (bw/1e6) << std::endl;
    _usrp->set_rx_bandwidth(bw,ch);
    std::clog << "[UHD] " << boost::format("Actual RX Bandwidth: %f MHz...") % (_usrp->get_rx_bandwidth(ch)/1e6) << std::endl << std::endl;

    // system setup.
    boost::this_thread::sleep(boost::posix_time::seconds(SETUP_TIME)); //allow for some setup time
}

bool
Guhd::check_locked_sensor(std::vector<std::string> sensor_names, const char* sensor_name, get_sensor_fn_t get_sensor_fn){
    if (std::find(sensor_names.begin(), sensor_names.end(), sensor_name) == sensor_names.end())
        return false;

    boost::system_time start = boost::get_system_time();
    boost::system_time first_lock_time;

    std::clog << boost::format("Waiting for \"%s\": ") % sensor_name;
    std::clog.flush();

    while (true) {
        if ((not first_lock_time.is_not_a_date_time()) and
                (boost::get_system_time() > (first_lock_time + boost::posix_time::seconds(SETUP_TIME)))) {
            std::clog << " locked." << std::endl;
            break;
        }
        if (get_sensor_fn(sensor_name).to_bool()){
            if (first_lock_time.is_not_a_date_time())
                first_lock_time = boost::get_system_time();
            std::clog << "+";
            std::clog.flush();
        } else {
            first_lock_time = boost::system_time();	//reset to 'not a date time'

            if (boost::get_system_time() > (start + boost::posix_time::seconds(SETUP_TIME))){
                std::clog << std::endl;
                throw std::runtime_error(str(boost::format("timed out waiting for consecutive locks on sensor \"%s\"") % sensor_name));
            } (std::clog << "_").flush();
        }
        boost::this_thread::sleep(boost::posix_time::milliseconds(100));
    }
    //std::clog << std::endl;
    return true;
}

