/*
 * BedFile.cpp
 *
 *  Created on: Nov 27, 2010
 *      Author: zhangc
 */
#include "BedFile.h"

#include <vector>
#include <string>
//#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include <assert.h>
#include <stdlib.h>

//#include <boost/algorithm/string.hpp>
//#include <boost/algorithm/string/join.hpp>

//extern "C" {
//jksrc header
//#include "common.h"
//#include "sqlNum.h"
//#include "utils.h"
//}

//my own header
#include "utils.h"


using namespace std;
//using namespace boost;


BedLine &
BedLine::operator=(const BedLine & bl)
{
	chrom = bl.chrom;
	chromStart = bl.chromStart;
	chromEnd = bl.chromEnd;
	name = bl.name;
	score = bl.score;
	strand = bl.strand;
	thickStart = bl.thickStart;
	thickEnd = bl.thickEnd;
	itemRgb = bl.itemRgb;
	blockCount = bl.blockCount;
	blockSizes = bl.blockSizes;
	blockStarts = bl.blockStarts;
	colNum = bl.colNum;
	return *this;
}

//parse a line to BedLine
BedLine::BedLine (const char *line)
{
	//vector <string> cols;
	//split (cols, line, is_any_of ("\t "));
	char cols[15][1024];

	//colNum = chopByWhite(const_cast <char *>(line), cols, 15);
	istringstream ss(line);
	string s;
	int ii = 0;
	for( ;ii<15; ii++ ) {
		if(getline(ss, s, '\t')) {
			sscanf ( s.c_str(), "%s", cols[ii] );
		}
		else break;
	}
	colNum = ii;

	if (colNum < 3)
	{
			cerr << "bed line must have at least 3 columns: " << line << endl;
			exit(1);
	}

	chrom = cols[0];
	chromStart = atoi(cols[1]);
	chromEnd = atoi(cols[2]) - 1;
	if (colNum > 3)
		name = cols[3];
	else
		return;

	if (colNum > 4)
		score = atof (cols[4]);
	else
		return;

	if (colNum > 5)
		strand = *cols[5];
	else
		return;

	if (colNum > 6)
	{
		assert (colNum > 7);
		thickStart = atoi (cols[6]);
		thickEnd = atoi (cols[7]) - 1;
	}
	else
		return;

	if (colNum > 8)
		itemRgb = cols[8];
	else
		return;

	if (colNum > 9)
	{
		blockCount = atoi (cols[9]);
		assert (blockCount > 0);
	}
	else
		return;

	if (blockCount > 1)
		assert (colNum > 11);

	if (colNum > 10)
	{
			assert (colNum > 11);
			//vector <string> v;
			string v[1024];
			ss.str( cols[10] );
			for( ii=0 ;ii<1024; ii++ ) {
				if(getline(ss, s, ',')) {
					v[ii] = s;
				}
				else break;
			}
			//size_t bCount = chopByChar (cols[10], ',', v, 1024);
			size_t bCount = ii;
			assert (bCount == blockCount);

			for (size_t i = 0; i < blockCount; i++)
				blockSizes.push_back (atoi(v[i].c_str() ));

			//split (v, cols[10], is_any_of (","));
			//blockSizes.resize(v.size());
			//transform (v.begin(), v.end (), blockSizes.begin(), string2int);
			
			//bCount = chopByChar (cols[11], ',', v, 1024);
			ss.str( cols[11]);
			for( ii=0 ;ii<1024; ii++ ) {
				if(getline(ss, s, ',')) {
					 v[ii] = s;
				}
				else break;
			}
			bCount = ii;
			assert (bCount == blockCount);

			for (size_t i = 0; i < blockCount; i++)
				blockStarts.push_back (atoi(v[i].c_str() ));

			//split (v, cols[11], is_any_of (","));
			//blockStarts.resize(v.size());
			//transform (v.begin(), v.end (), blockStarts.begin(), string2int);
			//assert (blockSizes.size() == blockStarts.size());
	}
}

//expand it to the 12 column format
void BedLine::BedLineToFull ()
{
	if (colNum <= 4)
	{
		stringstream ss;
		ss << chrom << ":" << chromStart << "-" << (chromEnd+1);
		name = ss.str();
	}

	if (colNum <= 5)
		score = 0;

	if (colNum <= 6)
		strand = '+';

	if (colNum <= 7)
		thickStart = chromStart;

	if (colNum <= 8)
		thickEnd = chromEnd;

	if (colNum <= 9)
		itemRgb = "0,0,0";

	if (colNum <= 10)
		blockCount = 1;

	if (colNum <= 11)
		blockSizes.push_back (chromEnd - chromStart + 1);

	if (colNum <= 12)
		blockStarts.push_back (0);
}

/*
BedLine::BedLine (const string &line)
{
	vector <string> cols;
	split (cols, line, is_any_of ("\t "));

	if (cols.size() < 3)
	{
			cerr << "bed line must have at least 3 columns: " << line << endl;
			exit(1);
	}

	colNum = cols.size();

	chrom = cols[0];
	chromStart = string2int(cols[1]);
	chromEnd = string2int(cols[2]) - 1;
	if (colNum > 3)
		name = cols[3];
	else
		return;

	if (colNum > 4)
		score = string2float (cols[4]);
	else
		return;

	if (colNum > 5)
		strand = cols[5].c_str()[0];
	else
		return;

	if (colNum > 6)
	{
		assert (colNum > 7);
		thickStart = string2int(cols[6]);
		thickEnd = string2int (cols[7]) - 1;
	}
	else
		return;

	if (colNum > 8)
		itemRgb = cols[8];
	else
		return;

	if (colNum > 9)
	{
		blockCount = string2int (cols[9]);
		assert (blockCount > 0);
	}
	else
		return;

	if (blockCount > 1)
		assert (colNum > 11);

	if (colNum > 10)
	{
			assert (colNum > 11);
			vector <string> v;

			split (v, cols[10], is_any_of (","));
			blockSizes.resize(v.size());
			transform (v.begin(), v.end (), blockSizes.begin(), string2int);

			split (v, cols[11], is_any_of (","));
			blockStarts.resize(v.size());
			transform (v.begin(), v.end (), blockStarts.begin(), string2int);
			assert (blockSizes.size() == blockStarts.size());
	}
}

*/
string BedLine::BedLineToString ()
{
	stringstream ss;

	ss << chrom << "\t" << chromStart << "\t" << chromEnd+1;

	if (colNum > 3)
		ss << "\t" << name;

	if (colNum > 4)
		ss << "\t" << score;

	if (colNum > 5)
		ss << "\t" << strand;

	if (colNum > 6)
	{
		ss << "\t" << thickStart;
		assert (colNum > 7);
		ss << "\t" << thickEnd + 1;
	}

	if (colNum > 8)
		ss << "\t" << itemRgb;

	if (colNum > 9)
		ss << "\t" << blockCount;

	if (colNum > 9 && blockCount > 1)
	{
		assert (colNum > 11);
	}

	if (colNum > 10)
	{
		ss << "\t" << join (blockSizes.begin(), blockSizes.end(), ",");
		assert (colNum > 11);
		ss << "\t" << join (blockStarts.begin(), blockStarts.end(), ",");
	}

	return ss.str();
}

size_t BedLine::GetColNum ()
{
	return colNum;
}

string BedLine::GetChrom ()
{
	assert (colNum >= 3);
	return chrom;
}

int64_t BedLine::GetChromStart ()
{
	assert (colNum >=3);
	return chromStart;
}

int64_t BedLine::GetChromEnd ()
{
	assert (colNum >=3);
	return chromEnd;
}

string BedLine::GetName ()
{
	assert (colNum >= 4);
	return name;
}

float BedLine::GetScore ()
{
	assert (colNum >= 5);
	return score;
}

char BedLine::GetStrand ()
{
	assert (colNum >= 6);
	return strand;
}

int64_t BedLine::GetThickStart ()
{
	assert (colNum >= 8);
	return thickStart;
}

int64_t BedLine::GetThickEnd ()
{
	assert (colNum >= 8);
	return thickEnd;
}

string BedLine::GetItemRgb ()
{
	assert (colNum >= 9);
	return itemRgb;
}

size_t BedLine::GetBlockCount ()
{
	assert (colNum >= 10);
	return blockCount;
}
const vector <int64_t> & BedLine::GetBlockSizes ()
{
	assert (colNum == 12 && blockSizes.size () == blockCount);
	return blockSizes;
}
const vector <int64_t> & BedLine::GetBlockStarts ()
{
	assert (colNum == 12 && blockSizes.size () == blockCount);
	return blockStarts;
}


BedFile::BedFile (const string &fn, bool v)
{
	filename = fn;
	verbose = v;
	in.open (filename.c_str());

	if (!in)
	{
		cerr << "cannot open file " << filename << " to read" << endl;
		exit (1);
	}
}

vector<BedLine>& BedFile::NextBedLineBlock (int maxGap = 0, bool separateStrand =false, size_t minBlockSize = 0)
{
	int64_t currBlockEnd = -1;
	bedLineBlock.clear ();

	//the max gap allowed between bed regions in a block
	assert (maxGap >= 0);
	streampos currPointer = in.tellg ();

	while (BedLine* pbl = NextBedLine ())
	{
		if (separateStrand)
			assert (pbl->GetColNum () >= 6);

		bool expand = false;

		if (currBlockEnd < 0 || bedLineBlock.size () < minBlockSize)
		{
			expand = true;
		}	
		else
		{
			size_t currSize = bedLineBlock.size ();
			BedLine lastBl = bedLineBlock[currSize-1];
			if (separateStrand)
			{
				if (lastBl.GetStrand () == pbl->GetStrand ()
					&& lastBl.GetChrom ().compare (pbl->GetChrom ()) == 0
					&& pbl->GetChromStart () - currBlockEnd -1 <= maxGap)
				{
					expand = true;
				}
			}
			else
			{
				if (lastBl.GetChrom ().compare (pbl->GetChrom ()) == 0
					&& pbl->GetChromStart () - currBlockEnd -1 <= maxGap)
				{
					expand = true;
				}
			
			}
		}

		if (expand)
		{
			bedLineBlock.push_back (*pbl);
			if (pbl->GetChromEnd () > currBlockEnd)
				currBlockEnd = pbl->GetChromEnd ();
			currPointer = in.tellg ();
		}
		else
		{
			in.seekg (currPointer);
			return bedLineBlock;
		}
		delete pbl;
	}
	return bedLineBlock;
}


BedLine* BedFile::NextBedLine ()
{
	char line[bufferSize];
	//string line = "";
	while (!in.eof())
	{
		in.getline (line, bufferSize);
		//getline (in, line);
		if (lineStartWith (line, "#") || lineStartWith (line, "track")
			|| lineStartWith (line, "browser") || emptyLine (line))
			continue;
		break;
	}

	if (emptyLine (line))
		return (BedLine *) NULL;
	else
		return new BedLine (line);
		//return bl;
}

void BedFile::AddBedLine (BedLine & bl)
{
	bedLines.push_back (bl);
}

void BedFile::ReadFile ()
{
	size_t iter = 0;
	while (BedLine *pbl=NextBedLine())
	{
		if (iter % 100000 == 0 && verbose)
			cout << iter << " ..." << endl;
		iter++;

		AddBedLine (*pbl);
		delete pbl;
	}
}



