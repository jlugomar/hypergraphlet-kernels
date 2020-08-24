#include "string.h"


string capitalize(const string &s)  {
    string temp;

// TO BE CODED!

    return temp;
}


vector<string> split(const string &s, char delim)  {
    vector<string> temp;

    if (s.size() > 0)  {
        string::size_type prev(0), curr(0);

        while ((curr = s.find(delim, curr)) != string::npos)  {
            temp.push_back(s.substr(prev, curr-prev));
            curr++;
            prev = curr;
        }

        temp.push_back(s.substr(prev));        
    }

    return temp;
}


string strip(const string &s)  {
    bool foundch(false);
    unsigned start_(0);

    while (start_<s.size() && !foundch)  {
        if (isspace(s.at(start_))) 
            start_++;
        else
            foundch = true;
    } 
    
    foundch = false;
    unsigned end_(s.size()-1);

    while (end_>start_ && !foundch)  {
        if (isspace(s.at(end_)))
            end_--;
        else
            foundch = true;
    }

    return s.substr(start_, end_-start_+1);
}

