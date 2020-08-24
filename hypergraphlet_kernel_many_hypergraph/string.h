/**
 * Auxiliary string functions which immitate some of the useful methods of
 * Ruby strings.
 *
 * See http://www.ruby-doc.org/docs/ProgrammingRuby/html/ref_c_string.html  
 *
 * Vladimir Vacic, vladimir@cs.ucr.edu
 * Computer Science and Engineering
 * University of California, Riverside
 *
 * Apr-2-2008
 */

#ifndef __STRING_H__
#define __STRING_H__

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <cctype>
#include <sstream>
#include <string>
#include <vector>
using namespace std;


/** */
string capitalize(const string&);

/** String tokenizer. Splits the argument string using the delimiter
 *  parameter. */
vector<string> split(const string&, char);

/** Removes white space chatacters from begining and end of the string. */
string strip(const string &s);

inline int to_i(const string &s)  {
    return atoi(s.c_str());
}

inline double to_f(const string &s)  {
    return atof(s.c_str());
}

inline string to_s(int t) {
    ostringstream s;
    s << t;
    return s.str();
}

#endif
