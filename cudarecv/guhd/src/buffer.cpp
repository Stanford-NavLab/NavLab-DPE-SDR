
#include <iostream>
#include <sstream>
#include <ctime>

#include <boost/filesystem.hpp>
#include "guhd.h"
#include "buffer.h"

//Guhd::Buffer::~Buffer () {}

/*Guhd::Buffer*const
Guhd::Buffer::init(){
    if (!(*_p_run) || _init) return this;
    _init = true;
    return this;
}*/

/*Guhd::Buffer*const
Guhd::Buffer::terminate(){
    if (!_init) return this;
    _init = false;
    return this;
}*/

void
Guhd::FileBuff::append (const std::vector<std::complex<MEM_TYPE> >& readbuf,const size_t& num_of_samp){
    UHD_ASSERT_THROW (*_p_run);
    _ofs.write((char*)&readbuf.front(),sizeof(std::complex<MEM_TYPE>)*num_of_samp);
}

//// FILE BUF

Guhd::FileBuff::FileBuff (
        const boost::program_options::variables_map& vm,
        const size_t& chan,
        const bool&   ltwo,
        const volatile bool*const p_run):
    _p_run  (p_run),
    _ltwo   (ltwo),
    _ip     (vm["ip"].as<std::vector<size_t>>()[chan]),
    _rate   (vm["rate"].as<double>()),
    _path   (vm["path"].as<boost::filesystem::path>()){}

int
Guhd::FileBuff::open(const std::time_t& t){
    //if (!*_p_run ) return -1;
    UHD_ASSERT_THROW (*_p_run);

    namespace fs = boost::filesystem;
    if (!fs::exists(_path)) fs::create_directories(_path);
    UHD_ASSERT_THROW(fs::is_directory(_path));

    auto tm = *std::gmtime(&t);
    std::ostringstream oss;

    oss << _path.string() << '/'
        << std::put_time(&tm, "%Y%m%d_%H%M%S")
        << "_usrp" << _ip
        << "_" << size_t(_rate/1e3) << "kHz"
        << ".dat";
    std::clog << "[UHD] Writing new file to: " << oss.str() << std::endl;

    if (_ofs.is_open()) _ofs.close();
    _ofs.open(oss.str(), std::ios::binary);

    const fs::path syml (_path.string() + "/usrp" + std::to_string(_ip) + ".dat");
    if (fs::exists(syml) && !fs::is_symlink(syml)) {
        std::cerr << "[UHD] No symbolic link created for latest file." << std::endl;
        return -1;
    }

    if (fs::is_symlink(syml)) fs::remove(syml);
    fs::create_symlink(oss.str(),syml);
    return 0;
}

int
Guhd::FileBuff::close(){
    if(_ofs.is_open()) {
        _ofs.close();
        std::clog << "[UHD] Closing file." << std::endl;
    } return 0;
}

