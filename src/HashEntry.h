/*
 * HashEntry.h
 *
 *  Created on: Jul 24, 2014
 *      Author: galaxy
 */
#include <string>
#include "ReadSequence.h"


#ifndef HASHENTRY_H_
#define HASHENTRY_H_

namespace HardSearch {

class HashEntry {
public:
	HashEntry(uint32_t hc);
	virtual ~HashEntry();

	std::string GetKey() { return this->_k; };
	uint32_t GetHashCode() { return this->_hc; };
	ReadSequence* GetValue() { return this->_val; };

	void SetValue(ReadSequence *new_value) { this->_k = new_value->readID; this->_val = new_value; };
	bool Equals(ReadSequence *other) { return this->_k == other->readID; };


private:

	/*
	 * @abstract			variables for hash entry
	 * @string		  _k	String storing key for the hash
	 * @uint32_t	  _hc	Int storing hashcode for the key
	 * @ReadSequence  _val	Pointer to the readsequence(value)
	 */
	std::string _k;
	uint32_t _hc;
	ReadSequence *_val;



};

} /* namespace HardSearch */

#endif /* HASHENTRY_H_ */
