#ifndef DOCONFIG_H
#define DOCONFIG_H


#include <vector>
#include <string>

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

#endif
