/*
 * main.cpp
 *
 *  Created on: Jul 29, 2014
 *      Author: galaxy
 */
#include <stdio.h>
#include <string>
#include <cstring>
#include <cstdlib>


int main(int argc, char *argv[])
{

	std::string cmd;
	std::string rs;
	std::string blast_path;
	std::string blast_ref;

	rs = std::string(argv[1], strlen(argv[1]));
	blast_path = std::string(argv[2], strlen(argv[2]));
	blast_ref = std::string(argv[3], strlen(argv[3]));

	cmd = 	"echo '" + rs	+ "' | " + blast_path +
			" -db " + blast_ref + " -num_threads 1 -max_hsps 1 " +
			"-max_target_seqs 1 -gapopen 1000 -gapextend 10 -outfmt 5";

	FILE *ls = popen(cmd.c_str(), "r");
	char* buf = new char[4096]();
	std::string stream_string;
	std::string xml_string;
	while (fgets(buf, sizeof(buf), ls) != 0) {
		stream_string = std::string(buf, strlen(buf));
		xml_string += stream_string;

	}
	pclose(ls);


	printf(xml_string.c_str());
	printf("\n");

	return(0);

}


