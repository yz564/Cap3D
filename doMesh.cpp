#include "doMesh.h"
#include "doConfig.h"
#include <assert.h>

template <typename T>
Mesh<T>::Mesh(int Ne, int Nn, int Et) : num_elem(Ne), num_node(Nn), elem_type(Et), elem_ptrs(Ne), nodes(Nn) {
	printf("[Mesh] constructor: \n elem_num=%d, node_num=%d \n", num_elem, num_node);
}

template <typename T>
T Triangle<T>::calArea(std::vector<Node<T>> * nodes_ptr){
	Node p1=(*nodes_ptr)[0];
	Node p2=(*nodes_ptr)[1];
	Node p3=(*nodes_ptr)[2];
	T a = std::sqrt(std::pow(p2.x - p1.x, 2) + std::pow(p2.y - p1.y, 2) + std::pow(p2.z - p1.z, 2));
    	T b = std::sqrt(std::pow(p3.x - p2.x, 2) + std::pow(p3.y - p2.y, 2) + std::pow(p3.z - p2.z, 2));
    	T c = std::sqrt(std::pow(p1.x - p3.x, 2) + std::pow(p1.y - p3.y, 2) + std::pow(p1.z - p3.z, 2));
    	T s = (a + b + c) / 2;  // Semi-perimeter
    	this->area = std::sqrt(s * (s - a) * (s - b) * (s - c));  // Heron's formula
    	return this->area;
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
	int num_blk = tmp[0]; int num_elem = tmp[1]; int num_node = tmp[2];
	get_line_strip_comments(readfile, str); //third line
	ReadLine_helper<int>(str, 1, tmp);
	int elem_type = tmp[0];
	Mesh<T> * mesh_ptr = new Mesh<T>(num_elem, num_node, elem_type);
	int count = 0;
	for (int i = 0; i < num_blk; ++i) {
		get_line_strip_comments(readfile, str);
		ReadLine_helper<int>(str, 3, tmp);
		int num_elem_blk = tmp[0];
		//int blk_id = tmp[1];
		int attrib_blk = tmp[2];
		if (elem_type==3){
			for (int j = 0; j < num_elem_blk; ++j) {
				Triangle<T> * elem_ptr=new Triangle<T>(j, attrib_blk);
				std::getline(readfile, str);
				ReadLine_helper<int>(str, elem_type, elem_ptr->node_ids);
				for (int k = 0; k < elem_type; ++k){elem_ptr->node_ids[k]-=1;} //meshfile index is 1-based
				mesh_ptr->elem_ptrs[count]=elem_ptr;
				count++;
			}
		}else{
			printf("Failed: The element type of %d nodes has not been implemented\n",elem_type);
			break;
		}
	}
	assert(count == num_elem);
	//start to read the node coordinates (three-dimensional)
	std::vector<T> tmp_t(3);
	for (int i = 0; i < num_node; ++i) {
		std::getline(readfile, str);
		ReadLine_helper<T>(str, 3, tmp_t);
		Node<T> tmp_node(i,tmp_t[0],tmp_t[1],tmp_t[2]);
		//mesh_ptr->nodes[i]=Node<T>{i,tmp_t[0],tmp_t[1],tmp_t[2]};
		mesh_ptr->nodes[i]=tmp_node;
	}
	//calculate area in each Element
	T total_area=0;
	std::vector<Node<T>> nodes_in_elem(elem_type);
	for (int i=0; i<num_elem; ++i){
		for (int j=0; j<elem_type; ++j){
			nodes_in_elem[j]=(mesh_ptr->nodes[mesh_ptr->elem_ptrs[i]->node_ids[j]]);
			(mesh_ptr->nodes[mesh_ptr->elem_ptrs[i]->node_ids[j]]).link_elem_ids.push_back(i);
		}
		total_area+=mesh_ptr->elem_ptrs[i]->calArea(&nodes_in_elem);
	}
	std::cout << "the total area of the mesh: " << total_area <<std::endl;
	return mesh_ptr;
}
//explicit instantiations for commonly used types
//otherwise it results in "undefined reference" errors in linker
template Mesh<float> * LoadMesh<float>(const std::string meshfile);
template Mesh<double> * LoadMesh<double>(const std::string meshfile);

template<typename T>
void Node<T>::print_info(std::ofstream & logfile){
	logfile<<"Node id:"<< node_id <<std::endl;
	logfile<<"Coordinates x: "<< x <<", y: "<< y <<", z: " << z <<std::endl;
	logfile<<"Elements that includes this node: ";
	for (size_t i=0; i<link_elem_ids.size(); ++i){
		logfile<< link_elem_ids[i] <<" ";
	}
	logfile<<std::endl<<std::endl;
}
template void Node<float>::print_info(std::ofstream & logfile);
template void Node<double>::print_info(std::ofstream & logfile);


