#ifndef DOCONFIG_H
#define DOCONFIG_H

#include <iostream>
#include <vector>
#include <string>
#include <fstream>

template <typename T>
class Config {
public:
	T relative_permitivity;
	int metal_num;
	std::vector<T> electrostatic_potential;
	int allowedthreads;
	const std::string mesh_file;
	
	Config(T eps_r, int Nm, std::vector<T> phi, int tds, const std::string mf);
	void print();
};

template <typename T>
Config<T> * ReadConfig(const std::string inputfile);

template<typename T>
void ReadLine_helper(const std::string & str, const int num, std::vector<T> & data);

std::istream& get_line_strip_comments(std::istream& stm, std::string& str);


#endif
