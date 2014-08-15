/*
 * SamParser.h
 *
 *  Created on: Jun 20, 2014
 *      Author: galaxy
 */
#include <string>
#include <iostream>
#include <fstream>
#include <exception>
#include <algorithm>
#include <cmath>
#include <map>
#include <pthread.h>

#include "SoftClipRead.h"
#include "BedUtils.h"
#include "sam.h"
#include "bam.h"
#include "Utils.h"
#include "MurmurHash3.h"
#include "HashMap.h"
#include "HashEntry.h"


#ifndef SAMPARSER_H_
#define SAMPARSER_H_



class SamParser {
public:

	/*
	 * @abstract	structure for bam fetch callback
	 */
	struct args
	{
		bam_header_t *bh;
		bam1_t *b;
		int quality;
		std::string direction;
		// TODO: Should be able to remove these variables...
		//std::vector<ReadSequence *> *unmapped_reads_ptr;
		//std::vector<std::string> *unmapped_reads_set;

		HardSearch::HashMap *unmapped_reads_ptr2;
	};

	struct parse_args
	{
		pthread_mutex_t *mut;
		ReadSequence *first_read;
		ReadSequence *second_read;
		std::vector<bool>* thread_status;
		std::vector<SoftClipRead *> *softclip_reads_ptr;
		int min_softclip_size;

		std::string reference_fp;
		int current_thread;
		int counter;

		//std::vector<struct parse_args *>* jobs;
	};

	struct job_list
	{
		pthread_mutex_t *mut;
		int unsigned completed_jobs;
		std::vector<struct parse_args *> jobs;
		std::vector<int> failed_jobs;

		uint64_t *left_softclip_reads;
		uint64_t *right_softclip_reads;
	};


	//std::vector<SoftClipRead> softclip_reads;
	std::vector<ReadSequence> discordant_reads;
	long unsigned total_softclip_reads;
	long unsigned left_softclip_reads;
	long unsigned right_softclip_reads;
	uint32_t mean;
	uint32_t deviation;


	std::map< std::string, std::vector<ReadSequence> > target_reads_ptr;
	std::vector<SoftClipRead *> softclip_reads_ptr;
	std::vector<ReadSequence *> deletion_reads_ptr;
	std::vector<ReadSequence *> unmapped_reads_ptr;
	std::vector<std::string> unmapped_reads_set;

	HardSearch::HashMap discordant_reads_ptr;
	HardSearch::HashMap target_reads_ptr2;
	HardSearch::HashMap unmapped_reads_ptr2;

	/*
	 * @abstract			Create the Sam Parser object for parsing BAM/SAM files
	 * 						using the SAMTools API
	 * @param  min_stdev	Float specifying the minimum standard deviation from mean
	 * 						to consider discordant
	 * @param  min_sc_size	Int (unsigned) specifying the minimum softclip size to
	 * 						consider for target region
	 * @param  min_rs_qual	Int specifying the minimum mapping quality of the read to
	 * 						consider for target region
	 *
	 * @discussion	Sam Parser class provides a wrapper for opening and parsing SAM files
	 * using the SAMTools API. In addition, the class contains a list of all the softclip
	 * and discordant reads used for further analysis downstream
	 */
	SamParser(	float min_stdev,
				uint32_t min_sc_size,
				uint32_t min_rs_qual);

	SamParser();
	virtual ~SamParser();

	/*
	 * @absract			Opens the sam/bam file and parses it for softclip and discordant reads
	 * @param  sam_fp	String specifying the file path of the BAM file
	 * @param  ref_fp	String specifying the file path of the reference genome
	 *
	 * @discussion	Opens up the BAM file and parses the reads using the SAMTools api. Reads with
	 * softclip or discordant reads are collected for use downstream.
	 */
	void Parse(std::string sam_fp, std::string ref_fp);

	void static ParseForSoftClips(	ReadSequence rs,
									pthread_mutex_t *mut,
									uint64_t &left_counter,
									uint64_t &right_counter,
									int min_softclip_size,
									std::string reference_fp,
									std::vector<SoftClipRead *> *softclip_reads_ptr);

	void ParseForDeletions(ReadSequence rs, pthread_mutex_t *mut);
	void ParseForUnmappedReads(ReadSequence rs);


	// set the path and open the sam file for reading
	// Parse the provided Sam file and store all the reads
	// containing soft clips
	void SetParameters(	float min_standard_deviation,
						uint32_t min_softclip_size,
						uint32_t min_mapping_quality);

	/*
	 * @abstract	Print out the softclip read statistics
	 *
	 * @discussion	Outputs the number of total softclips, left softclips,
	 * and right softclips found in the BAM file to the stdout
	 */
	void PrintNumberOfSoftClips();

	/*
	 * @abstract	Print out the mean and stdev for the alignments
	 *
	 * @discussion	Outputs the mean and standard deviation of the inserte
	 * sizes of alignments in the BAM file. The statistics only take into account
	 * reads that have insert sizes less than 1000 bps currently (there is an
	 * issue where outlying discordant read pairs are throwing off the statistics)
	 */
	void PrintAlignmentStatistics();


	/*
	 * @abstract	Fix left and right softclip reads
	 * @param  rs	SoftClipRead to fix
	 *
	 * @discussion	Softclip reads that do not get called properly by the aligner
	 * must be fixed before generating an8cefoprstud calling target regions
	 */
	void FixSoftClip(SoftClipRead rs);

	/*
	 * @abstract		Write out all end clips in BED format
	 * @param  out_fn	String containing the filename to write out to
	 *
	 * @discussion		Used to output all fixed end clips for
	 */
	bool OutputSoftClipReads(std::string out_fn);

	//
	// warpper methods for grabbing alignment information from BAM file
	//
	/*
	 * @abstract	Grab the read id from the alignment
	 * @param  b	Pointer to the structure of the alignment
	 * @return		String containing the read id
	 */
	static std::string GetReadID(const bam1_t *b);

	/*
	 * @abstract	Grab the bitflag field from the alignment
	 * @param  b	Pointer to the structure of the alignment
	 * @return		Int containing the bitflag
	 */
	static uint32_t GetBitflag(const bam1_t *b);

	/*
	 * @abstract	Grab the chromosome field from the alignment
	 * @param  b	Pointer to the structure of the alignment
	 * @param  bh	Pointer to the header of the BAM file
	 * @return		String containing the chromosome
	 */
	static std::string GetChromosome(const bam1_t *b, bam_header_t *bh);

	/*
	 * @abstract	Grab the start coordinate from the alignment
	 * @param  b	Pointer to the structure of the alignment
	 * @return		Int containing the coordinate of the left most
	 * 				aligning base
	 */
	static int32_t GetStartCoordinate(const bam1_t *b);

	/*
	 * @abstract	Grab the cigar field from the alignment
	 * @param  b	Pointer to the structure of the alignment
	 * @return		String containing the cigar field
	 */
	static std::string GetCigar(const bam1_t *b);

	/*
	 * @abstract	Grab the quality score field from the alignment
	 * @param  b	Pointer to the structure of the alignment
	 * @return		Int containing the quality score
	 */
	static uint32_t GetQualityScore(const bam1_t *b);

	/*
	 * @abstract	Grab the read sequence field from the alignment
	 * @param  b	Pointer to the structure of the alignment
	 * @return		String containing the read sequence
	 */
	static std::string GetReadSequence(const bam1_t *b);

	/*
	 * @abstract	Grab the insert size of the alignment from the
	 * 				template length field
	 * @param b		Pointer to the structure of the alignment
	 * @return		Int containing the template size
	 */
	static uint32_t GetInsertSize(const bam1_t *b);

	/*
	 * @abstract				Set the alignment statistics (std, mean) insert
	 * 							size for the SAM file
	 *
	 * @discussion	calculate the standard deviation and mean of the insert size
	 * of all the reads in the SAM file
	 */
	void _GetAlignmentStatistics();


	/*
	 * @abstract		Get reads corresponding to the target regions
	 * @param  tid		Int specifying target id from header
	 * @param  beg		Int specifying the lower bound to search
	 * @param  end		Int specifying the upper bound to search
	 * @param  direct	String specifying the direction of region to search
	 *
	 * @discussion		Wrapper for SAMTools api function to retrieve reads
	 * for given target region
	 */
	void RetrieveReadsForRegion(int tid,
								int beg,
								int end,
								std::string direct);

	void _GetMateReadID(std::string sam_file_path, HardSearch::HashMap reads);

	std::string GetSamFilePath() { return this->_sam_file_path; };

private:

	uint32_t _min_softclip_size;
	uint32_t _min_mapping_quality;
	float _min_standard_deviation;
	std::string _reference_file_path;
	std::string _sam_file_path;

	float _upper_limit;
	float _lower_limit;

	/*
	 * @abstract				Get all the mates for target region reads
	 * @param  sam_file_path	String containing the path to the alignment file
	 *
	 * @discussion	Grab and set all of the mates for Softclip and Discordant Reads
	 */
	void _GetMateReadID(std::string sam_file_path);

	static ReadSequence _GetReadFromAlignment(const bam1_t *b, bam_header_t *bh);
	static int _UnmappedMatesCallBack(const bam1_t *b, void *data);
	/*
	 * @abstract	 parse the soft clip to determine which type of read it is
	 */
	void _InitDefaultValues();
	void _ParseReads();
	static void *_Parse(void *arguments);
};

#endif /* SAMPARSER_H_ */
