
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include "flowmgr.h"
#include "auxil.h"

dsp::FlowMgr::Registrar::~Registrar(){}

dsp::FlowMgr::FlowMgr(){
    _regis.reserve(5);
    _regis.push_back(new dsp::FlowMgr::DPE);
}

dsp::FlowMgr::~FlowMgr(){
    for (auto const& id: _flowlist)
        if (id.first) delete id.first;
    for (auto const& ptr:_regis) delete ptr;
}

size_t
dsp::FlowMgr::EmergencyStop()const{
    std::cout << "[FlorMgr] Emergency Stop. All flow stopping." << std::endl;
    size_t actual = 0;
    for (auto& pair: _flowlist){
        if (!pair.first) continue;
        if (pair.first->Stop() == 0) actual ++;
    } return actual;
}

dsp::FlowMgr::FlowMgr(const FlowMgr& rhs):_alias(rhs._alias){
    _flowlist.reserve(rhs._flowlist.size());
    for (auto const& pair: rhs._flowlist)
        createFlow(pair.second);
}//hopefully will never need this.

int
dsp::FlowMgr::createFlow(const std::string& type, const std::string& alias){
    FlowIdent fid (nullptr,"");

    for(auto const& ptr:_regis){
        if (!ptr->match(type)) continue;
        fid.first  = ptr->createFlow();
        fid.second = ptr->get_name();
    }

    if (!fid.first){
        std::cerr << "[FlowMgr] Flow type \"" << type << "\" does not exist." << std::endl;
        return -1;
    }

    _flowlist.push_back(fid);
    const size_t flowN = _flowlist.size() - 1;
    if (!alias.empty()) addAlias(alias,flowN);
    std::clog << "[FlowMgr] Flow #"<< flowN << " (" << fid.second << ") created." << std::endl;
    return 0;
}

int
dsp::FlowMgr::destroyFlow(const size_t& idx){
    if (idx >= _flowlist.size())
        throw std::out_of_range("[FlowMgr] destroyFlow: idx exceeded max range.");

    if (!_flowlist[idx].first){
        std::cerr<< "[FlowMgr] Flow #" << idx << " previously destroyed." << std::endl;
        return -1;
    }
    //_flowlist[idx].first -> Stop();
    delete _flowlist[idx].first;
    _flowlist[idx].first = nullptr;
    std::clog << "[FlowMgr] Flow #" << idx <<  " (" << _flowlist[idx].second <<") succesfully deleted." << std::endl;
    return 0;
}

int
dsp::FlowMgr::startFlow(const size_t& idx)const{
    try{
        int status = -1;
        if ((status = getFlowPtr(idx)->Start())!=0) return status;
        std::clog << "[FlowMgr] Flow #" << idx <<  " started." << std::endl;
    } catch (const std::out_of_range& e) {
        std::cerr << e.what() << std::endl;
    } return 0;
}

int
dsp::FlowMgr::stopFlow(const size_t& idx)const{
    try{
        int status = -1;
        if ((status = getFlowPtr(idx)->Stop()) != 0) return status;
        std::clog << "[FlowMgr] Flow #" << idx <<  " being stopped." << std::endl;
    } catch (const std::out_of_range& e) {
        std::cerr << e.what() << std::endl;
    } return 0;
}

int
dsp::FlowMgr::loadFlow(const std::string& key,const char* fname)const{
    size_t idx = getFlowIdx (key);
    return (idx == NPOS)?-1:loadFlow(idx,fname);
}

int
dsp::FlowMgr::loadFlow(const size_t& idx,const char* fname)const{
    int status = -1;
    try {
        status = getFlowPtr(idx)->LoadFlow(fname);
        std::clog << "[FlowMgr] Flow #" << idx <<  " being loaded." << std::endl;
    } catch (const std::out_of_range& e) {
        std::cerr << e.what() << std::endl;
    } return status;
}

size_t
dsp::FlowMgr::getFlowIdx(const std::string& key)const{
    size_t idx = NPOS;

    bool allnum = true;
    for (auto const& c: key){
        if (!isalnum(c)){
            std::cerr << "[FlowMgr] Alias \"" << key << "\" contains invalid character."<<std::endl;
            return NPOS;
        } allnum = (allnum && isdigit(c));
    }

    if (!key.empty() && allnum){
        size_t idx = size_t(atoi(key.c_str()));
        if (idx >= _flowlist.size()) return NPOS;
        return idx;
    }

    try {
        idx = _alias.at (key);
    } catch (const std::out_of_range& e){
        std::cerr <<"[FlowMgr] Alias \"" << key <<"\" does not exist." <<std::endl;
        return NPOS;
    } return idx;
}

dsp::Flow*
dsp::FlowMgr::getFlowPtr(const size_t& idx)const{
    if (idx >= _flowlist.size())
        throw std::out_of_range("[FlowMgr] getFlow: idx exceeded max range.");
    dsp::Flow* ptr = _flowlist[idx].first;
    if (!ptr) throw std::logic_error("[FlowMgr] Flow constructor not properly invoked.");
    return _flowlist[idx].first;
}

void
dsp::FlowMgr::addAlias (const std::string& key,const size_t& value){
    if (value >= _flowlist.size()){
        std::cerr << "[FlowMgr] addAlias: idx exceeded max range." << std::endl;
        return;
    } else if (key.size() > 80 || key.empty()){
        std::cerr << "[FlowMgr] Maximum length for alias is 80 characters." << std::endl;
        return;
    }

    for (auto const& c:key)
        if (!isalnum(c)) {
            std::cerr << "[FlowMgr] Alias \"" << key <<"\" invalid; must be alphanumerals." << std::endl;
            return;
        }

    if (!isalpha(key[0]))
        std::cerr << "[FlowMgr] Alias \"" << key <<"\" invalid; first char must be alphabet." << std::endl;
    else if (_alias.insert(std::pair<std::string,size_t>(key,value)).second)
        std::clog << "[FlowMgr] Alias \"" << key <<"\" created for Flow #" << value << std::endl;
    else
        std::cerr << "[FlowMgr] Alias \"" << key <<"\" already existed for Flow #" << _alias[key] << '.' << std::endl;
}

void
dsp::FlowMgr::listAlias()const{
    if (_alias.empty()){
        std::cout << "No active alias." << std::endl;
        return;
    }

    for (auto const & entry: _alias){
        const size_t i = entry.second;
        if (!_flowlist[i].first) continue;
        std::cout << std::setw (15) <<std::left << ("\"" + entry.first + "\" -> ");
        std::cout << std::setw(3) << i << " : ("
                  << std::setw(15) << std::left << _flowlist[i].second
                  << ") " << _flowlist[i].first << std::endl;
    }
}

void
dsp::FlowMgr::listFlow()const{
    if (_flowlist.empty()){
        std::cout << "No active flow." << std::endl;
        return;
    }
    for(size_t i = 0 ; i < _flowlist.size() ;i++){
        if (!_flowlist[i].first) continue;
        std::cout << std::setw(3) << i << " : ("
                  << std::setw(15) << std::left << _flowlist[i].second
                  << ") " << _flowlist[i].first << std::endl;
    }
}

int
dsp::FlowMgr::setParam (const std::string& key, const std::string& modName, const std::string& paramName, const std::string& expr)const{
    size_t idx = getFlowIdx (key);
    return (idx == NPOS)?-1:setParam(idx,modName,paramName,expr);
}

int
dsp::FlowMgr::setParam (const size_t& idx, const std::string& modName, const std::string& paramName, const std::string& expr)const{

    std::clog << "[FlowMgr] Param type = " ;
    if (/*expr.size() == 3 && expr[0] == '\'' && expr[1] != '\'' && expr[2] == '\''*/ expr[0] == '\\' && expr[1] == 'x'){
        std::clog << "char" << std::endl;
        char val = std::stoi(expr.substr(2),nullptr,16);
        return getFlowPtr(idx) -> SetModParam (modName,paramName,val);
    } else if (expr[0] == '"' && expr[expr.size() - 1] == '"'){
        const std::string val = expr.substr(1,expr.size() -2);
        if (val.find('"') == std::string::npos) {
            std::clog << "char*" << std::endl;
            return getFlowPtr(idx) -> SetModParam(modName,paramName,val.c_str());
        } std::clog << "UNDEFINED" << std::endl; return -1;
    } else if (expr == "true") {
        std::clog << "bool" << std::endl;
        return getFlowPtr(idx) -> SetModParam(modName,paramName,true);
    } else if (expr == "false") {
        std::clog << "bool" << std::endl;
        return getFlowPtr(idx) -> SetModParam(modName,paramName,false);
    } else if (expr.find_first_not_of("0123456789.+-ef") != std::string::npos){
        std::clog << "bool" << std::endl;
        return -1;
    }

    const size_t sign = expr.find_first_of ("+-");
    if (sign != 0 && sign != std::string::npos) {
        std::clog << "UNDEFINED" << std::endl;
        return -1;
    }

    std::istringstream iss (expr);

    size_t decimal = expr.find_first_of(".ef");
    if (decimal != std::string::npos){
        if (decimal == (expr.size() - 1) || expr.find('.',decimal + 1) == std::string::npos){ // float
            float val;
            iss >> val;
            std::clog << "float" << std::endl;
            return getFlowPtr(idx) -> SetModParam(modName,paramName,val);
        } std::clog << "UNDEFINED" << std::endl; return -1;
    }

    int val;
    iss >> val;
    std::clog << "int" << std::endl;
    return getFlowPtr(idx) -> SetModParam(modName,paramName,val);
}

// "list" and "monitor"

int
dsp::FlowMgr::listOutput (const std::string& key, const std::string& modName,const std::string& portName)const{
    size_t idx = getFlowIdx (key);
    return (idx == NPOS)?-1:listOutput(idx,modName,portName);
}

int
dsp::FlowMgr::listOutput(const size_t& idx, const std::string& modName,const std::string& portName)const{
    dsp::Port* ptr = nullptr;
    try {
        if (getFlowPtr(idx) -> GetOutput(modName,portName,&ptr) != 0 || !ptr)
            throw std::out_of_range("[FlowMgr] Cannot retrieve designated port.");
    } catch (const std::out_of_range& e) {
        std::cerr << e.what() << std::endl;
        return -1;
    }

    size_t size = ptr->VectorLength;

    if (!size){
        std::cout << "Empty." << std::endl;
        return 0;
    }

    switch(ptr->Datatype){
        /*case FLOAT_t:
            size *= sizeof(float);
            std::cout << "Float: [";
            break;
        case DOUBLE_t:
            size *= sizeof(double);
            std::cout <<" Double [";
            break;
        case INT_t:
            size *= sizeof(int);
            std::cout << "Int [";
            break;*/
        case CHAR_t:
            size *= sizeof(char);
            std::cout <<" Char [";
            break;
        case BOOL_t:
            size *= sizeof(bool);
            std::cout << "Bool [";
            break;
        default:
            std::cerr << "[FlowMgr] Data type unhandled." <<std::endl;
            return -1;
    }

    void* data = malloc(size);
    memcpy(data,ptr->Data,size);

    if (ptr->Datatype == BOOL_t){
        bool* bdata = (bool*) data;
        for(size_t i = 0 ; i < ptr->VectorLength ; i++){
            if (bdata[i]) std::cout << i << ',';
        } std::cout << "\b]" << std::endl;
    } else if (ptr->Datatype == CHAR_t){
        char* cdata = (char*) data;
        for(size_t i = 0 ; i < ptr->VectorLength ; i++)
            std::cout << (int)cdata[i] << ',';
        std::cout << "\b]" <<std::endl;
    }   free (data);
    return 0;
}

void
dsp::FlowMgr::flowType()const{
    for(auto const& ptr:_regis)
        std::cout << std::setw(14) << std::left << ptr->get_name() << ' '<< ptr->get_desc() << std::endl;
}
