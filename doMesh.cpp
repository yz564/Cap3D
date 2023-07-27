#include "doMesh.h"
#include "doConfig.h"
#include <assert.h>

template <typename T>
Mesh<T>::Mesh(int Ne, int Nn, int Et, int Na) : num_elem(Ne), num_node(Nn), elem_type(Et), num_attrib(Na), elem_ptrs(Ne), nodes(Nn) {
	printf("[Mesh] constructor: Elem_num=%d, node_num=%d \n", num_elem, num_node);
	logfile<<"[Mesh] constructor: Elem_num = " <<num_elem<<", node_num = " <<num_node<< std::endl;
}

template <typename T>
T Triangle<T>::calArea(){
	Node<T>* p1=this->node_ptrs[0];
	Node<T>* p2=this->node_ptrs[1];
	Node<T>* p3=this->node_ptrs[2];
	T a = std::sqrt(std::pow(p2->x - p1->x, 2) + std::pow(p2->y - p1->y, 2) + std::pow(p2->z - p1->z, 2));
    	T b = std::sqrt(std::pow(p3->x - p2->x, 2) + std::pow(p3->y - p2->y, 2) + std::pow(p3->z - p2->z, 2));
    	T c = std::sqrt(std::pow(p1->x - p3->x, 2) + std::pow(p1->y - p3->y, 2) + std::pow(p1->z - p3->z, 2));
    	T s = (a + b + c) / 2.0;  // Semi-perimeter
    	this->area = std::sqrt(s * (s - a) * (s - b) * (s - c));  // Heron's formula
    	return this->area;
}
// need to modify for non-parallel cases
template <typename T>
T Quadrilateral<T>::calArea(){
	Node<T>* p1=this->node_ptrs[0];
	Node<T>* p2=this->node_ptrs[1];
	Node<T>* p3=this->node_ptrs[2];
	//Node<T>* p4=this->node_ptrs[3];
	Node<T> ab= Node<T>(-1,p2->x-p1->x,p2->y-p1->y,p2->z-p1->z);
	Node<T> ac= Node<T>(-1,p3->x-p1->x,p3->y-p1->y,p3->z-p1->z);
	Node<T> N= Node<T>(-1, ab.y*ac.z-ab.z*ac.y, ab.z*ac.x-ab.x*ac.z, ab.x*ac.y-ab.y*ac.x);
	this->area = std::sqrt(N.x * N.x + N.y * N.y + N.z * N.z)/1.0;
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
	Mesh<T> * mesh_ptr = new Mesh<T>(num_elem, num_node, elem_type, num_blk);
	int count = 0;
	for (int i = 0; i < num_blk; ++i) {
		get_line_strip_comments(readfile, str);
		ReadLine_helper<int>(str, 3, tmp);
		int num_elem_blk = tmp[0];
		//int blk_id = tmp[1];
		int attrib_blk = tmp[2]-1;
		if (elem_type==3){
			for (int j = 0; j < num_elem_blk; ++j) {
				Triangle<T> * elem_ptr=new Triangle<T>(j, attrib_blk);
				std::getline(readfile, str);
				ReadLine_helper<int>(str, elem_type, elem_ptr->node_ids);
				for (int k = 0; k < elem_type; ++k){elem_ptr->node_ids[k]-=1;} //meshfile index is 1-based
				mesh_ptr->elem_ptrs[count]=elem_ptr;
				count++;
			}
		}else if(elem_type==4){
			for (int j = 0; j < num_elem_blk; ++j) {
				Quadrilateral<T> * elem_ptr=new Quadrilateral<T>(j, attrib_blk);
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
	
	std::vector<T> total_area(num_blk,0);
	int degree = std::stoi(std::getenv("QUADRATURE_DEGREE"));
        printf("Quadrature degree is %d \n", degree);
	//calculate area of each Element
	//std::vector<Node<T>> nodes_in_elem(elem_type);
	for (int i=0; i<num_elem; ++i){
		Element<T>* temp_elem = mesh_ptr->elem_ptrs[i];
		for (int j=0; j<elem_type; ++j){
			temp_elem->node_ptrs[j]=&(mesh_ptr->nodes[temp_elem->node_ids[j]]);
			//mesh_ptr->elem_ptrs[i]->node_ptrs[j]=&mesh_ptr->nodes[mesh_ptr->elem_ptrs[i]->node_ids[j]];
			//nodes_in_elem[j]=(mesh_ptr->nodes[mesh_ptr->elem_ptrs[i]->node_ids[j]]);
			temp_elem->node_ptrs[j]->link_elem_ids.push_back(i);
			//(mesh_ptr->nodes[mesh_ptr->elem_ptrs[i]->node_ids[j]]).link_elem_ids.push_back(i);
		}
		total_area[temp_elem->attrib]+=temp_elem->calArea();
		temp_elem->setQuadrature(1);
		temp_elem->center=temp_elem->qpoints[0]; //calculate the center point;
		temp_elem->qpoints.pop_back();
		assert(temp_elem->qpoints.empty());
		temp_elem->setQuadrature(degree);
	}
	for (int i=0; i<num_blk; ++i){
		std::cout << "The area of metal objective " << i+1 << ": " << total_area[i] << "m^2 \n";
		logfile << "The area of metal objective " << i+1 << ": " << total_area[i] << "m^2 \n";
	}
	return mesh_ptr;
}
//explicit instantiations for commonly used types
//otherwise it results in "undefined reference" errors in linker
template Mesh<float> * LoadMesh<float>(const std::string meshfile);
template Mesh<double> * LoadMesh<double>(const std::string meshfile);

template<typename T>
void Node<T>::print_info(){
	logfile<<"Node id:"<< node_id <<std::endl;
	logfile<<"Coordinates x: "<< x <<", y: "<< y <<", z: " << z <<std::endl;
	logfile<<"Elements that includes this node: ";
	for (size_t i=0; i<link_elem_ids.size(); ++i){
		logfile<< link_elem_ids[i] <<" ";
	}
	logfile<<std::endl<<std::endl;
}
template void Node<float>::print_info();
template void Node<double>::print_info();


