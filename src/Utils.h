/*
 * Utils.h
 *
 *  Created on: Jul 17, 2014
 *      Author: galaxy
 */
#include <string>
#include <sstream>
#include <iomanip>


#ifndef UTILS_H_
#define UTILS_H_

namespace HardSearch {

template<typename T>
std::string NumToString(T number, int unsigned precision = 0)
{
	std::string sResult;
	std::ostringstream tempStream;
	tempStream << std::fixed << std::setprecision(precision) << number;
	sResult = tempStream.str();
	return sResult;
}

class Utils {
public:
	Utils();
	virtual ~Utils();
};

} /* namespace HardSearch */

#endif /* UTILS_H_ */
