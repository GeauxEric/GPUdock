//Header for mpidecoys.C

#include <iostream>
#include <new>
#include <fstream>
#include <string>
#include "complex.h"
#include <limits>
#include <vector>

using namespace std;

void Decoy(Complex *, std::string, std::string);

std::ifstream& GotoLine(std::ifstream&, int );
