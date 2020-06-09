#ifndef _STREAMBLOCK_H
#define _STREAMBLOCK_H

#include <boost/format.hpp>
#include <boost/thread.hpp>
#include <boost/lexical_cast.hpp>

#include <vector>
#include <string>
#include <complex>
#include <fstream>
#include <cassert>

#include <uhd/types/tune_request.hpp>
#include <uhd/utils/thread_priority.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <boost/program_options.hpp>

#ifdef SC8
    #define STREAMBLOCK_MEM_TYPE char
    #define STREAMBLOCK_OTW_FORMAT "sc8"
    #define STREAMBLOCK_CPU_FORMAT "sc8"
#else
    #define STREAMBLOCK_MEM_TYPE short
    #define STREAMBLOCK_OTW_FORMAT "sc16"
    #define STREAMBLOCK_CPU_FORMAT "sc16"
#endif


class dsp::streamblock {
    using get_sensor_fn_t = boost::function<uhd::sensor_value_t (const std::string&)>;

    volatile bool               _run,_overflow_message;
    const std::string           _clock,_sync;
    const bool                  _mimo;
    const double                _rate;

    uhd::usrp::multi_usrp::sptr _usrp;
    uhd::rx_streamer::sptr      _rx_stream;
    std::vector<size_t>         _chan;
    std::vector<bool>           _b_ltwo;
    size_t                      _spb,_runN;

    std::vector<std::vector<std::complex<MEM_TYPE> > > _buffs;
    std::vector<std::complex<MEM_TYPE> *>              _buff_ptrs;
    std::vector<dsp::streamblock::CircBuff<std::vector<std::complex<MEM_TYPE> > > > _circ_buffs;
    
    

    std::string get_fname   (const std::time_t&,const size_t&)const;
    void        set_clock   ();
    void        set_rf      (const size_t&,const bool = false);

    static bool
    check_locked_sensor     (std::vector<std::string> sensor_names, const char* sensor_name, get_sensor_fn_t get_sensor_fn);

    class Buffer;
    class FileBuff;
    class SockBuff;
    std::vector<FileBuff*>  _p_fbuff;
    //std::vector<SockBuff*>  _p_sock;

public:
    dsp::streamblock(
/*        const std::vector<size_t>& ip,
        const std::vector<size_t>& ltwo,
        const std::string& clk,
        const std::string& sync,*/
        const boost::program_options::variables_map& vm,
        const bool&   mimo/*,
        const double& rate,
        const boost::filesystem::path& path *///= boost::filesystem::current_path()
    );

    ~dsp::streamblock();

    int Start();
    int Update();
    int Stop();
};

#endif
