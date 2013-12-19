//Header for mpidecoys.C

#include <iostream>
#include <new>
#include <fstream>
#include <string>
#include <limits>
#include <vector>
#include "complex.h"

using namespace std;

void Decoy(Complex *, std::string, std::string);

std::ifstream& GotoLine(std::ifstream&, int );
