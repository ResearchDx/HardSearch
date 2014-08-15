/*
 * RapidXMLWrapper.h
 *
 *  Created on: Jul 17, 2014
 *      Author: galaxy
 */
#include "XmlParser/rapidxml.hpp"
#include "XmlParser/rapidxml_utils.hpp"
#include "XmlParser/rapidxml_print.hpp"
#include "ReadSequence.h"
#include "Utils.h"
#include <string>
#include <map>

#ifndef RAPIDXMLWRAPPER_H_
#define RAPIDXMLWRAPPER_H_

namespace HardSearch {

class RapidXMLWrapper {
public:

	rapidxml::xml_document<>* doc;
	std::map<int, rapidxml::xml_node<>*> *doc_nodes;


	RapidXMLWrapper();
	virtual ~RapidXMLWrapper();

	// Given the name of a node in string format, this will
	// append the node to the xml document and add a pointer
	// to the node in the doc_nodes map for access
	rapidxml::xml_node<>* AddNode(std::string node_name);

	/*
	 * @abstract		Add an attribute to the document
	 * @param  name	String containing the name of the attribute
	 * @param  value	String containing the value of the attribute
	 * @return		Pointer to the xml_attribute

	 * @discussion	Allocates an xml_attribute to the document. Returns a
					pointer to the attribute with the provided name and value
	 */


	rapidxml::xml_attribute<>* AddAttribute(std::string attribute_name, std::string value);

	// Return true of the node exists in the map of nodes in
	// the current document, otherwise return false
	bool NodeExists(std::string node_name);

	rapidxml::xml_node<>* AddReadNode(ReadSequence &rs);
	void AddReadAlignmentToNode(AlignmentNode &ba, rapidxml::xml_node<>* node);
	void AddReadDeletionToNode();

	bool OutputXml(std::string fn);
	static void SetNodeValue(rapidxml::xml_node<>* &node, std::string value);

	// Get the node with given name
//	rapidxml::xml_node<>* GetNode(std::string node_name);


private:
	int _node_count;

	// Append the version and encoding to the start of the xml string
	void _InitializeXMLDocument();




};

} /* namespace HardSearch */

#endif /* RAPIDXMLWRAPPER_H_ */
