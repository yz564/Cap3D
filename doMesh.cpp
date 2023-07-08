#include "doMesh.h"
#include "doConfig.h"
#include <assert.h>

template <typename T>
Mesh<T>::Mesh(int Ne, int Nn, int Et) : num_elem(Ne), num_node(Nn), elem_type(Et), elems(Ne), attribs(Ne), nodes(Nn) {
	for (int i=0; i< num_elem; ++i){
		elems[i].resize(elem_type);
	}
	for (int i=0; i< num_node; ++i){
		nodes[i].resize(3);
	}
	printf("[Mesh] constructor: \n elem_num=%d, node_num=%d \n", num_elem, num_node);
}

/* 
this function should change if the mesh template changes.
current mesh template are in the ./mesh/xxx.elx
typename T can be float or double
*/
template<typename T>
Mesh<T> * LoadMesh(const std::string meshfile) {
	std::ifstream readfile(meshfile);
	std::string str;
	//first line is a comment
	get_line_strip_comments(readfile, str); //second line has  num_blk  num_elem  num_node
	std::vector<int> tmp(3);
	ReadLine_helper<int>(str, 3, tmp);
	int num_blk = tmp[0];
	int num_elem = tmp[1];
	int num_node = tmp[2];
	get_line_strip_comments(readfile, str); //third line
	ReadLine_helper<int>(str, 1, tmp);
	int elem_type = tmp[0];
	Mesh<T> * ans = new Mesh<T>(num_elem, num_node, elem_type);
	int count = 0;
	for (int i = 0; i < num_blk; ++i) {
		get_line_strip_comments(readfile, str);
		ReadLine_helper<int>(str, 3, tmp);
		int num_elem_blk = tmp[0];
		//int blk_id = tmp[1];
		int attrib_blk = tmp[2];
		for (int j = 0; j < num_elem_blk; ++j) {
			std::getline(readfile, str);
			ReadLine_helper<int>(str, elem_type, ans->elems[count]);
			ans->attribs[count]=attrib_blk;
			count++;
		}
	}
	assert(count == num_elem);
	//start to read the node coordinates (three-dimensional)
	for (int i = 0; i < num_node; ++i) {
		std::getline(readfile, str);
		ReadLine_helper<T>(str, 3, ans->nodes[i]);
	}
	return ans;
}
//explicit instantiations for commonly used types
//otherwise it results in "undefined reference" errors in linker
template Mesh<float> * LoadMesh<float>(const std::string meshfile);
template Mesh<double> * LoadMesh<double>(const std::string meshfile);



