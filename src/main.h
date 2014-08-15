/*	S&W Alignment Post SoftSearch
*/
#ifndef LIB_INCLUDES
#define LIB_INCLUDES
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <stdio.h>
#include <cstdlib>
#include <unistd.h>
#include <dirent.h>
#include <pthread.h>
#include <exception>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/thread/thread.hpp>
#include <ctime>
#include <sys/wait.h>
#include <sys/types.h>


#include "XmlParser/rapidxml.hpp"
#include "XmlParser/rapidxml_utils.hpp"
#include "XmlParser/rapidxml_print.hpp"

#include "ReadSequence.h"
#include "AlignmentNode.h"
#include "AlignmentStructs.h"
#include "BreakPointNode.h"
#include "BedUtils.h"
#include "SamParser.h"
#include "RapidXMLWrapper.h"
#include "Utils.h"
#include "MurmurHash3.h"
#include "HashMap.h"

#endif

#ifndef MAIN_H
#define MAIN_H

namespace FS = boost::filesystem;

int main(int argc, char *argv[]);

template<typename T>
std::string NumToString(T number)
{
	std::string sResult;
	std::ostringstream tempStream;
	tempStream << std::fixed << std::setprecision(0) << number;
	sResult = tempStream.str();
	return sResult;
}

template<typename T>
T GetParam(boost::program_options::variable_value val)
{
	try
	{
		return val.as<T>();
	}
	catch (boost::bad_any_cast &e)
	{
		printf("\n\tParameter Value Error\n");
		throw boost::bad_any_cast();
	}
}


struct arg_struct {
	std::string read_sequence;
	bool mate_mode;
	std::vector<ReadSequence>* read_list;
	pthread_mutex_t *mut;
	std::vector<bool>* thread_status;
	int current_thread;
	int counter;

};

struct job_queue{
	pthread_mutex_t *mut;
	uint32_t completed_jobs;
	std::vector<struct arg_struct *>* jobs;
	std::vector<int> failed_jobs;
	std::string path_to_blasthelper;
	std::string path_to_blast_reference;
	std::string path_to_blast;
};

std::string CheckFilePath(	boost::program_options::variable_value val, FS::path cwd,
							std::string param, bool blast_database =  false);

// @abstract	Print the help text to the screen
void DisplayHelp();

void _AddProgramParameters(boost::program_options::options_description &desc);
void _ParseProgramParameters(boost::program_options::variables_map vm);

int SizeOfParameterList(char** argv);

// Locates softclips in the file /output/directory/hardsearch.direction.bed
// Merges softclips spanning the same event using sortBed command
// Softclips are organized by direction
// Outputs the result to /output/directory/hardsearch.fixed.bed
void MergeSoftClips(std::string output_directory);

// Construct string from merged softclip region
std::string _MergeSoftClipRegion(	std::vector<std::string> softclips,
									std::vector<std::string> softclip_data,
									std::string _tid,
									std::string direction);

// Get the genomic coordinates of regions with greater than min_soft_clipped_count
// soft clip reads mapped to it, these will be the regions of interest when
// searching for additional unmapped softclip reads
std::vector<std::string> GetTargetRegions(std::string input_file, int min_soft_clipped_count);
void GetAdditionalTargetRegions(std::ifstream &target_region_file, std::vector<std::string> &target_regions);

int GetSoftClippedDirection(std::string soft_clipped_read);

/*
 * @abstract		Worker function that blasts reads in queue
 * @param  args		Void Pointer to job_queue struct that contains
 *					all the information for each job and queue
 *
 * @discussion		Iterates through the job queue and forks a call to
 *					BlastHelper. Output from blast helper is piped back
 *					in xml format. Output is parsed and the alignment
 *					information is set for the read
 */
void *_BlastReads(void *args);


/*
 * @abstract		Parse a XML node and store in ReadSequence
 * @param  node		Pointer to xml_node with read alignment information
 * @param  r		ReadSequence to set data for
 * @param mm		Bool specifying if alignment for read mate
 *
 * @dicussion		Parse alignment data in XML format from Blast and
 * 					set the corresponding ReadSequence with the results
 */
void ParseNode(rapidxml::xml_node<> *node, ReadSequence& r, bool mm);

// This function will split a line by a delimiter and
// return the result in a vector
static std::vector<std::string>SplitByDelimiter(char *input, char *delimiter);

std::string _ConstructFastaFileMates(ReadSequence rs);

static std::string ConstructFastaRead(std::string readName, std::string readSequence);

// Get the directory of a file path
std::string GetDirectory(char *filepath);

void ParseBlastXML(std::string inputFile, std::vector<ReadSequence> &readList, ReadSequence &rs, bool mateMode);

bool CompareChr(ReadSequence read1, ReadSequence read2);

// Return the path of the executable
std::string get_selfpath();

// Write out any reads that have an unmapped mate in the target region
bool _RetrieveUnmappedMateReads(std::string path_to_sam,
								std::string output_directory,
								std::vector<std::string> target_regions,
								int unsigned mapping_quality,
								SamParser *sp);

#endif
