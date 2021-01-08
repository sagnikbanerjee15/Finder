/*
 * utils.cpp
 *
 *  Created on: Nov 27, 2010
 *      Author: zhangc
 */

#include <string>
#include <sstream>

//c header
#include <string.h>
#include <stdlib.h>

//my own header
#include "utils.h"


using namespace std;

bool lineStartWith (const string &line, const string &pattern)
{
		return strncmp (line.c_str(), pattern.c_str(), pattern.size()) == 0 ? true : false;
}

bool emptyLine (const string &line)
{
		return line.find_first_not_of("\t ") == string::npos ? true : false;
}


int string2int(const string &str)
{
	return atoi(str.c_str());
}

float string2float (const string &str)
{
	return atof(str.c_str());
}


