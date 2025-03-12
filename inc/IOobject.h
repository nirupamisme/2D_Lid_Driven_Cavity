#ifndef __IOobject_H
#define __IOobject_H

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <map>
#include "IOobject.h"
#include "UserTypes.h"

using namespace std;

// Class for writing the data
class WriteData {
private:
    string filename;
    ofstream file;

public:
    WriteData(const string&);

    // Write function
    void write(array2D&);

    ~WriteData();
};

// Class for reading the input data
class ReadData {
private:
    string filename;
    ifstream file;

public:
    ReadData(const string&);

    // Read function
    void read(map<string, string>&);

    ~ReadData();
};

#endif      // __IOobject_H
