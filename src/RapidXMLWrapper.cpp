/*
 * RapidXMLWrapper.cpp
 *
 *  Created on: Jul 17, 2014
 *      Author: galaxy
 */

#include "RapidXMLWrapper.h"

namespace HardSearch {

RapidXMLWrapper::RapidXMLWrapper() {
	this->doc = new rapidxml::xml_document<>;
	this->doc_nodes = new std::map<int, rapidxml::xml_node<>*>;
	this->_node_count = 0;
	_InitializeXMLDocument();
}


void RapidXMLWrapper::_InitializeXMLDocument()
{
	rapidxml::xml_node<>* declaration = (*this->doc).allocate_node(rapidxml::node_declaration);
	rapidxml::xml_attribute<>* version = (*this->doc).allocate_attribute("version", "1.0");
	rapidxml::xml_attribute<>* encoding = (*this->doc).allocate_attribute("encoding", "utf-8");

	declaration->append_attribute(version);
	declaration->append_attribute(encoding);

	(*this->doc).append_node(declaration);
}

rapidxml::xml_node<>* RapidXMLWrapper::AddNode(std::string node_name)
{

	char *stored_node_name = (*this->doc).allocate_string(node_name.c_str());
	rapidxml::xml_node<>* node = (*this->doc).allocate_node(rapidxml::node_element, stored_node_name);
	this->_node_count += 1;
	(*this->doc_nodes)[this->_node_count] = node;
	return (*this->doc_nodes)[this->_node_count];
}

rapidxml::xml_attribute<>* RapidXMLWrapper::AddAttribute(std::string attribute_name, std::string value)
{
	char *stored_attribute_name = (*this->doc).allocate_string(attribute_name.c_str());
	char *stored_attribute_value = (*this->doc).allocate_string(value.c_str());
	rapidxml::xml_attribute<>* attribute = (*this->doc).allocate_attribute(stored_attribute_name, stored_attribute_value);
	return attribute;
}

bool RapidXMLWrapper::NodeExists(std::string node_name)
{
	for(int index = 1; index < this->_node_count; index++)
	{
		if((*this->doc_nodes)[index]->name() == node_name.c_str())
			return true;
	}
	return false;
}

rapidxml::xml_node<>* RapidXMLWrapper::AddReadNode(ReadSequence &rs)
{
	rapidxml::xml_node<>* read_node = this->AddNode("Read");

	rapidxml::xml_node<>* read_id = AddNode("Read_ID");
	rapidxml::xml_node<>* read_bitflag = this->AddNode("Flag");
	rapidxml::xml_node<>* read_ref_name = this->AddNode("Reference_Name");
	rapidxml::xml_node<>* read_start_pos = this->AddNode("Start_Coordinate");
	rapidxml::xml_node<>* read_mapq = this->AddNode("Mapping_Quality");
	rapidxml::xml_node<>* read_cigar = this->AddNode("Cigar");
	rapidxml::xml_node<>* read_insert_size = this->AddNode("Insert_Size");
	rapidxml::xml_node<>* read_sequence = this->AddNode("Read_Sequence");
	rapidxml::xml_node<>* read_quality_string = this->AddNode("Quality");

	RapidXMLWrapper::SetNodeValue(read_id, rs.readID);
	RapidXMLWrapper::SetNodeValue(read_bitflag, HardSearch::NumToString<uint32_t>(rs.bitflag));
	RapidXMLWrapper::SetNodeValue(read_ref_name, rs.chromosome);
	RapidXMLWrapper::SetNodeValue(read_start_pos, HardSearch::NumToString<int32_t>(rs.start_pos));
	RapidXMLWrapper::SetNodeValue(read_mapq, HardSearch::NumToString<uint32_t>(rs.qual_score));
	RapidXMLWrapper::SetNodeValue(read_insert_size, HardSearch::NumToString<int32_t>(rs.insert_size));
	RapidXMLWrapper::SetNodeValue(read_cigar, rs.cigar);
	RapidXMLWrapper::SetNodeValue(read_sequence, rs.sequence);
	// read quality string

	read_node->append_node(read_id);
	read_node->append_node(read_bitflag);
	read_node->append_node(read_ref_name);
	read_node->append_node(read_start_pos);
	read_node->append_node(read_mapq);
	read_node->append_node(read_cigar);
	read_node->append_node(read_insert_size);
	read_node->append_node(read_sequence);
	read_node->append_node(read_quality_string);

	return read_node;
}

void RapidXMLWrapper::AddReadAlignmentToNode(AlignmentNode &ba, rapidxml::xml_node<>* node)
{
	rapidxml::xml_node<>* read_alignment = this->AddNode("Alignments");

	rapidxml::xml_node<>* node_softclip_start_2_1 = this->AddNode("Start_Coordinate_1");
	rapidxml::xml_node<>* node_softclip_orientation_2_1 = this->AddNode("Orientation_1");
	rapidxml::xml_node<>* node_softclip_match_2_1 = this->AddNode("Match_Coordinates_1");
	rapidxml::xml_node<>* node_softclip_bitscore_2_1 = this->AddNode("Bitscore_1");
	rapidxml::xml_node<>* node_softclip_length_2_1 = this->AddNode("Length_1");
	rapidxml::xml_node<>* node_softclip_flag_2_1 = this->AddNode("Flag_1");

	rapidxml::xml_node<>* node_softclip_start_2_2 = this->AddNode("Start_Coordinate_2");
	rapidxml::xml_node<>* node_softclip_orientation_2_2 = this->AddNode("Orientation_2");
	rapidxml::xml_node<>* node_softclip_match_2_2 = this->AddNode("Match_Coordinates_2");
	rapidxml::xml_node<>* node_softclip_bitscore_2_2 = this->AddNode("Bitscore_2");
	rapidxml::xml_node<>* node_softclip_length_2_2 = this->AddNode("Length_2");
	rapidxml::xml_node<>* node_softclip_flag_2_2 = this->AddNode("Flag_2");


	std::string start_position_2_1 = ba.chromosome +":";
	start_position_2_1 += HardSearch::NumToString<uint64_t>(ba.aligned_start) + "-";
	start_position_2_1 += HardSearch::NumToString<uint64_t>(ba.aligned_end);

	std::string match_position_2_1 = HardSearch::NumToString<uint32_t>(ba.match_start) + "-";
	match_position_2_1 += HardSearch::NumToString<uint32_t>(ba.match_end);
	HardSearch::RapidXMLWrapper::SetNodeValue(node_softclip_start_2_1, start_position_2_1);
	HardSearch::RapidXMLWrapper::SetNodeValue(node_softclip_orientation_2_1, ba.Orientation());
	HardSearch::RapidXMLWrapper::SetNodeValue(node_softclip_match_2_1, match_position_2_1);
	HardSearch::RapidXMLWrapper::SetNodeValue(node_softclip_bitscore_2_1, HardSearch::NumToString<float>(ba.bit_score));
	HardSearch::RapidXMLWrapper::SetNodeValue(node_softclip_length_2_1, HardSearch::NumToString<uint32_t>(ba.length));
	HardSearch::RapidXMLWrapper::SetNodeValue(node_softclip_flag_2_1, HardSearch::NumToString<uint32_t>(ba.bitflag));

	if(ba.mate_alignment != NULL)
	{
		std::string start_position_2_2 = ba.mate_alignment->chromosome +":";
		start_position_2_2 += HardSearch::NumToString<uint64_t>(ba.mate_alignment->aligned_start) + "-";
		start_position_2_2 += HardSearch::NumToString<uint64_t>(ba.mate_alignment->aligned_end);

		std::string match_position_2_2 = HardSearch::NumToString<uint32_t>(ba.mate_alignment->match_start) + "-";
		match_position_2_2 += HardSearch::NumToString<uint32_t>(ba.mate_alignment->match_end);
		HardSearch::RapidXMLWrapper::SetNodeValue(node_softclip_start_2_2, start_position_2_2);
		HardSearch::RapidXMLWrapper::SetNodeValue(node_softclip_orientation_2_2, ba.mate_alignment->Orientation());
		HardSearch::RapidXMLWrapper::SetNodeValue(node_softclip_match_2_2, match_position_2_2);
		HardSearch::RapidXMLWrapper::SetNodeValue(node_softclip_bitscore_2_2, HardSearch::NumToString<float>(ba.mate_alignment->bit_score));
		HardSearch::RapidXMLWrapper::SetNodeValue(node_softclip_length_2_2, HardSearch::NumToString<uint32_t>(ba.mate_alignment->length));
		HardSearch::RapidXMLWrapper::SetNodeValue(node_softclip_flag_2_2, HardSearch::NumToString<uint32_t>(ba.mate_alignment->bitflag));

	}

	node->append_node(read_alignment);
	read_alignment->append_node(node_softclip_start_2_1);
	read_alignment->append_node(node_softclip_orientation_2_1);
	read_alignment->append_node(node_softclip_match_2_1);
	read_alignment->append_node(node_softclip_bitscore_2_1);
	read_alignment->append_node(node_softclip_length_2_1);
	read_alignment->append_node(node_softclip_flag_2_1);

	read_alignment->append_node(node_softclip_start_2_2);
	read_alignment->append_node(node_softclip_orientation_2_2);
	read_alignment->append_node(node_softclip_match_2_2);
	read_alignment->append_node(node_softclip_bitscore_2_2);
	read_alignment->append_node(node_softclip_length_2_2);
	read_alignment->append_node(node_softclip_flag_2_2);
}

bool RapidXMLWrapper::OutputXml(std::string fn)
{
	std::ofstream of;
	of.open(fn.c_str());
	of << (*this->doc);
	of.flush();
	of.close();
	return true;
}



void RapidXMLWrapper::SetNodeValue(rapidxml::xml_node<>* &node, std::string value) {
	std::string *ptr_Value;
	ptr_Value = new std::string(value);
	node->value(ptr_Value->c_str());
}

RapidXMLWrapper::~RapidXMLWrapper() {
	// delete XML document pointer to free up memory
	delete this->doc;
	delete this->doc_nodes;
}

} /* namespace HardSearch */
