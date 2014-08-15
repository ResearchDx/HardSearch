/*
 * HashEntry.cpp
 *
 *  Created on: Jul 24, 2014
 *      Author: galaxy
 */

#include "HashEntry.h"

namespace HardSearch {

HashEntry::HashEntry( uint32_t hc) {
	this->_hc = hc;
	this->_val = NULL;
}


HashEntry::~HashEntry() {
	// TODO Auto-generated destructor stub
}

} /* namespace HardSearch */
