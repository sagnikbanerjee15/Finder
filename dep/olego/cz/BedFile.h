/*
 * BedFile.h
 *
 *  Created on: Nov 27, 2010
 *      Author: zhangc
 */

#ifndef BEDFILE_H_
#define BEDFILE_H_

#include <vector>
#include <string>
#include <fstream>
#include <stdint.h>

using namespace std;
//using namespace boost;

class BedLine
{
	static const size_t bufferSize = 1024;

	string chrom;
	int64_t chromStart;
	int64_t chromEnd;
	string name;
	float score;
	char strand;
	int64_t thickStart;
	int64_t thickEnd;
	string itemRgb;
	size_t blockCount;
	vector<int64_t> blockSizes;
	vector<int64_t> blockStarts;

	size_t colNum;

public:
	BedLine () : colNum (0) {}
	//BedLine (const string &);
	BedLine (const char *);
	BedLine &operator=(const BedLine &);
	string BedLineToString ();
	void BedLineToFull ();

	size_t GetColNum ();

	string GetChrom ();
	int64_t GetChromStart ();
	int64_t GetChromEnd ();
	string GetName ();
	float GetScore ();
	char GetStrand ();
	int64_t GetThickStart ();
	int64_t GetThickEnd ();
	string GetItemRgb ();
	size_t GetBlockCount ();
	const vector <int64_t> & GetBlockSizes ();
	const vector <int64_t> & GetBlockStarts ();
};




class BedFile
{
	static const size_t bufferSize = 1024;

	string filename;
	vector <BedLine> bedLines; //To store all lines of the file
	vector <BedLine> bedLineBlock; //To store the current bed line block

	ifstream in;

	bool verbose;


public:
	BedFile () : filename(""), verbose(false) {}
	BedFile (const string &fn, bool v=false);
	BedLine* NextBedLine ();

	vector<BedLine>& NextBedLineBlock (int maxGap, bool separateStrand, size_t minBlockSize);

	void AddBedLine (BedLine &);

	//BedLine parseBedLine (const string &);


	void ReadFile ();
};



#endif /* BEDFILE_H_ */
