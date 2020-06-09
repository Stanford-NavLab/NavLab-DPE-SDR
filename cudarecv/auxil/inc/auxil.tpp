#ifndef INC__AUXIL_TPP_
#define INC__AUXIL_TPP_

#include <string>
#include <map>
#include <iostream>

//namespace auxil {

    template <typename T>
    T EnumParser<T>::ParseSomeEnum(const std::string &value) 
    { 
        typename std::map <std::string, T>::const_iterator iValue = enumMap.find(value);
        if (iValue  == enumMap.end()) {
        	//std::cerr << "[Auxil] ParseSomeEnum fail, got: " << value
            //          << std::endl;
		}
		else {
			//std::cout << "[Auxil] ParseSomeEnum found: " << value << std::endl;
    	}
    	return iValue->second;
    }
    
//}

#endif
