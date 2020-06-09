#ifndef _STREAMBLOCK_BUF_H
#define _STREAMBLOCK_BUF_H

#include <boost/filesystem.hpp>
//#include <sys/socket.h>
//#include <netinet/in.h>
#include <boost/program_options.hpp>
#include <boost/asio.hpp>
#include <boost/circular_buffer.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>
#include <boost/thread/thread.hpp>
#include <boost/call_traits.hpp>
#include <boost/bind.hpp>

#include <boost/timer/timer.hpp> // for auto_cpu_timer

#include <vector>
#include <complex>

// The base FileBuff class
class dsp::streamblock::FileBuff {
    const volatile bool* const      _p_run;
    const bool                      _ltwo;
    const size_t                    _ip,_rate;
    const boost::filesystem::path   _path;
    std::ofstream                   _ofs;

    public:
    FileBuff (
            const boost::program_options::variables_map& vm,
            const size_t& chan,
            const bool&   ltwo,
            const volatile bool*const p_run);
    ~FileBuff  () {}
    void  append  (const std::vector<std::complex<MEM_TYPE> >&,const size_t&) ;
    int open  (const std::time_t& = 0);
    int close ();
};


// Template class for bounded circular buffer
// Refer to: https://www.boost.org/doc/libs/1_61_0/doc/html/circular_buffer/examples.html
template <class T>
class dsp::streamblock::CircBuff {
public:

	typedef boost::circular_buffer<T> _container_type;
	typedef typename container_type::size_type _size_type;
	typedef typename container_type::value_type _value_type;
	typedef typename boost::call_traits<value_type>::param_type _param_type;

	explicit CircBuff(size_type capacity) : _m_unread(0), _m_container(capacity) {}

	void push_front(typename boost::call_traits<_value_type>::param_type item);

	void pop_back(value_type* pItem);
	
private:
	CircBuff(const CircBuff&);              // Disabled copy constructor.
	CircBuff& operator = (const CircBuff&); // Disabled assign operator.

	bool is_not_empty() const { return _m_unread > 0; }
	bool is_not_full() const { return _m_unread < _m_container.capacity(); }

	size_type _m_unread;
	container_type _m_container;
	boost::mutex _m_mutex;
	boost::condition _m_not_empty;
	boost::condition _m_not_full;
	
}; //


#endif
