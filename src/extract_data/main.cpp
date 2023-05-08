#include <fstream>
#include <iostream>

int readEOP(const std::string& filename, double**& polarMotion);

int main(int argc, char const* argv[]) {
  std::string filename = "data/earth-orientation-paramters/eop_2023.txt";
  double** polarMotion;
  int numlines = readEOP(filename, polarMotion);
  if (numlines == -1) {
    std::cout << "Error reading file" << std::endl;
    return -1;
  }
  std::cout << "Number of lines: " << numlines << std::endl;
  // print data as an array to copy in C
  std::cout << "double polarMotion[" << numlines << "][6] = {" << std::endl;
  for (int i = 0; i < numlines; i++) {
    std::cout << "  {";
    for (int j = 0; j < 7; j++) {
      std::cout << polarMotion[i][j];
      if (j < 6) {
        std::cout << ", ";
      }
    }
    std::cout << "}";
    if (i < numlines - 1) {
      std::cout << ",";
    }
    std::cout << std::endl;
  }
  return 0;
}

// read earth-orbiting-parameters data from data/eop/
int readEOP(const std::string& filename, double**& polarMotion) {
  // for each line in the file, read the first 6 columns (yymmdd), the columns 19-27 for the x-coordinate of polar motion, the columns 38-46 for the y-coordinate of polar motion, the columns 59-68 for the UT1-UTC.
  std::ifstream file(filename);

  if (!file.is_open()) {
    std::cout << "Error opening file" << std::endl;
    return -1;
  }

  // count number of lines in file
  int n_lines = 0;
  std::string line;
  while (std::getline(file, line)) {
    n_lines++;
  }

  polarMotion = new double*[n_lines];
  for (int i = 0; i < n_lines; i++) {
    polarMotion[i] = new double[7];  // year, month, day, x, y, UT1-UTC
  }
  // reset file pointer to beginning of file
  file.clear();
  file.seekg(0, std::ios::beg);
  for (int i = 0; i < n_lines; i++) {
    std::getline(file, line);
    std::string year = line.substr(0, 2);
    std::string month = line.substr(2, 2);
    std::string day = line.substr(4, 2);
    std::string mjd = line.substr(7, 8);
    std::string x = line.substr(18, 9);
    std::string y = line.substr(37, 9);
    std::string UT1_UTC = line.substr(58, 9);

    polarMotion[i][0] = std::stod(year);
    polarMotion[i][1] = std::stod(month);
    polarMotion[i][2] = std::stod(day);
    polarMotion[i][3] = std::stod(mjd);
    polarMotion[i][4] = std::stod(x);
    polarMotion[i][5] = std::stod(y);
    polarMotion[i][6] = std::stod(UT1_UTC);
  }
  return n_lines;
}