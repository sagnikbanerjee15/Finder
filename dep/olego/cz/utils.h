/*
 * utils.h
 *
 *  Created on: Nov 27, 2010
 *      Author: zhangc
 */

#ifndef UTILS_H_
#define UTILS_H_

//c++ header
#include <string>

//c header


using namespace std;

bool lineStartWith (const string & /*line*/, const string & /*pattern*/);
bool emptyLine (const string &/*line*/);
int string2int(const string &/*str*/);
float string2float (const string &/*str*/);

//template <class T> string join(const T &/*begin*/, const T &/*end*/, const string &/*t*/);
//template <class T> 
//string join(const T &, const T &, const string &);

template <class T>
string join(const T &begin, const T &end, const string &t)
{
	stringstream ss;
	for (T it=begin; it!=end; it++)
	{
		if (!ss.str().empty())
			ss << t;
		ss << *it;
	}
	return ss.str();
}


#endif /* UTILS_H_ */
