#include "IOobject.h"

using namespace std;

// WriteData Class Functions
WriteData::WriteData(const string& f): filename(f), file(filename) {}

// Write function
void WriteData::write(array2D& u) {
    for (int i = 0; i < u.size(); i++) {
        for (int j = 0; j < u[i].size(); j++) {
            if (j == 0)
                file << u[j][i];
            else
                file << "," << u[j][i];
        }
        file << endl;
    }
}

WriteData::~WriteData() {
    file.close();
}

// ReadData Class Functions
ReadData::ReadData(const string& f): filename(f), file(filename) {}

// Read function
void ReadData::read(map<string, string>& m) {
    string line;
    int pos;
    while (getline(file, line)) {
        pos = line.find('=');
        if (pos != string::npos) {
            string key = line.substr(0, pos);
            string value = line.substr(pos + 1);
            m[key] = value;
        }
    }
}

ReadData::~ReadData() {
    file.close();
}
