#ifndef DOMESH_H
#define DOMESH_H

#include <array>
#include <sstream>
#include <iostream>
#include <vector>
#include <cstring>
#include <string>
#include <fstream>
#include <cmath>

template<typename T>
class Mesh {
public:
	int num_elem;
	int num_node;
	int elem_type;
	std::vector<std::vector<int>> elems;
	std::vector<int> attribs;
	std::vector<std::vector<T>> nodes;

	Mesh(int Ne, int Nn, int Et);
	
};



template<typename T>
Mesh<T> * LoadMesh(const std::string meshfile);

template<typename T>
void ReadLine_helper(const std::string & str, const int num, std::vector<T> & data);

std::istream& get_line_strip_comments(std::istream& stm, std::string& str);

#endif
