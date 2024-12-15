// Fucntions to load and procees pgm images.
// Author: Emmanuel A. Larralde Ortiz | ealarralde@gmail.com

#ifndef PGM_H
#define PGM_H

#include <string>
#include <vector>

typedef std::vector< std::vector<unsigned char> > Img;

Img read_pgm(std::string& fpath);
void write_pgm(std::string& fpath, Img& img);

class Pgm{
private:
    Img img;
public:
    Pgm(std::string& fpath);
    void print();
    std::vector<float> get_transformed();
};

void print(Pgm& pgm);

#endif