/*
 * HashMap.h
 *
 *  Created on: Jul 24, 2014
 *      Author: galaxy
 */

#include "MurmurHash3.h"
#include "HashEntry.h"
#include "ReadSequence.h"


#ifndef HASHMAP_H_
#define HASHMAP_H_

namespace HardSearch {

class HashMap {
public:

	uint32_t size;
	uint32_t item_count;

	HashMap();
	virtual ~HashMap();

	void Add(ReadSequence *rs);
	ReadSequence* Get(ReadSequence *rs);
	std::vector< std::vector< HashEntry* >* > GetMap() { return this->_map; };

private:
	std::vector< std::vector< HashEntry* >* > _map;
	int _GetHash(std::string rid);

};

} /* namespace HardSearch */

#endif /* HASHMAP_H_ */
