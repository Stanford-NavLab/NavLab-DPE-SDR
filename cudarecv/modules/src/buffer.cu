
#include <iostream>
#include <sstream>
#include <ctime>

#include <boost/filesystem.hpp>
#include "streamblock.h"
#include "buffer.h"


// *** FileBuff member functions ***

//dsp::Buffer::~Buffer () {}

/*dsp::Buffer*const
dsp::Buffer::init(){
    if (!(*_p_run) || _init) return this;
    _init = true;
    return this;
}*/

/*dsp::Buffer*const
dsp::Buffer::terminate(){
    if (!_init) return this;
    _init = false;
    return this;
}*/

void dsp::FileBuff::append (const std::vector<std::complex<MEM_TYPE> >& readbuf,const size_t& num_of_samp){
    UHD_ASSERT_THROW (*_p_run);
    _ofs.write((char*)&readbuf.front(),sizeof(std::complex<MEM_TYPE>)*num_of_samp);
}

//// FILE BUF

dsp::FileBuff::FileBuff (
        const boost::program_options::variables_map& vm,
        const size_t& chan,
        const bool&   ltwo,
        const volatile bool*const p_run):
    _p_run  (p_run),
    _ltwo   (ltwo),
    _ip     (vm["ip"].as<std::vector<size_t>>()[chan]),
    _rate   (vm["rate"].as<double>()),
    _path   (vm["path"].as<boost::filesystem::path>()){}

int dsp::FileBuff::open(const std::time_t& t){
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

    /*
    const fs::path syml (_path.string() + "/usrp" + std::to_string(_ip) + ".dat");
    if (fs::exists(syml) && !fs::is_symlink(syml)) {
        std::cerr << "[UHD] No symbolic link created for latest file." << std::endl;
        return -1;
    }
    if (fs::is_symlink(syml)) fs::remove(syml);
    fs::create_symlink(oss.str(),syml);
    */
    return 0;
}

int dsp::FileBuff::close(){
    if(_ofs.is_open()) {
        _ofs.close();
        std::clog << "[UHD] Closing file." << std::endl;
    } return 0;
}




// *** CircBuff Member Functions ***
// push_front() is called by the producer thread in order to insert a new item into the buffer.
// The method locks the mutex and waits until there is a space for the new item.
void dsp::CircBuff::push_front(typename boost::call_traits<value_type>::param_type item)
{ // `param_type` represents the "best" way to pass a parameter of type `value_type` to a method.
	
	boost::mutex::scoped_lock lock(m_mutex);
	m_not_full.wait(lock, boost::bind(&bounded_buffer<value_type>::is_not_full, this));
	m_container.push_front(item);
	++m_unread;
	lock.unlock();
	m_not_empty.notify_one();
}

// pop_back() method does NOT remove the item, but the item is left in the circular_buffer
// which then replaces it with a new one (inserted by a producer) when the circular_buffer is full.
// This technique is more effective than removing the item by explicitly calling the
// circular_buffer::pop_back() method of the circular_buffer, since assignment is more efficient than destruction
void dsp::CircBuff::pop_back(value_type* pItem) {
	boost::mutex::scoped_lock lock(m_mutex);
	m_not_empty.wait(lock, boost::bind(&bounded_buffer<value_type>::is_not_empty, this));
	*pItem = m_container[--m_unread];
	lock.unlock();
	m_not_full.notify_one();
}
