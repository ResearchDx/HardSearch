/*
 * HashMap.cpp
 *
 *  Created on: Jul 24, 2014
 *      Author: galaxy
 */

#include "HashMap.h"

namespace HardSearch {

HashMap::HashMap() {
	unsigned int i;

	for(i = 0; i < 2048; i++)
	{
		this->_map.push_back( new std::vector< HashEntry* > );
	}

	this->size = 2048;
	this->item_count = 0;

}

void HashMap::Add(ReadSequence *rs)
{
	uint32_t rid_hash = _GetHash((*rs).readID) % this->size;

	HashEntry *h = new HashEntry(rid_hash);
	h->SetValue(rs);

	if(this->Get(rs) == NULL)
	{
		this->_map[rid_hash]->push_back(h);
		this->item_count++;
	}
}

ReadSequence* HashMap::Get(ReadSequence *rs)
{
	uint32_t rid_hash = _GetHash((*rs).readID) % this->size;

	if(this->_map[rid_hash]->size() == 0)
		return NULL;

	std::vector< HashEntry* >::iterator h_iter;
	std::vector< HashEntry* > *h;

	h = this->_map[rid_hash];

	for(h_iter = h->begin(); h_iter != h->end(); h_iter++)
	{
		if((*h_iter)->Equals(rs))
			return (*h_iter)->GetValue();
	}

	return NULL;
}



int HashMap::_GetHash(std::string read_id)
{
	uint32_t hash_value = 0;
	MurmurHash3_x86_32(read_id.c_str(), strlen(read_id.c_str()), 0, &hash_value);

	return hash_value;
}


HashMap::~HashMap() {
	unsigned int i;
	for(i = 0; i < this->size; i++)
	{
		delete this->_map[i];
	}

}

} /* namespace HardSearch */
