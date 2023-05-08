#include "../include/misc.h"

double to2Pi(double angle) {
  return angle + 2 * PI * ((angle < 0) ? 1 : 0);
}

int numLines(std::string filename) {
  std::ifstream file(filename);
  int n = 0;
  std::string line;
  while (std::getline(file, line)) n++;
  // reset file and pointer to the beginning
  file.clear();
  file.seekg(0, std::ios::beg);
  return n;
}