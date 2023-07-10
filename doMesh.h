#ifndef DOMESH_H
#define DOMESH_H

#include <vector>
#include <string>
#include <fstream>
#include <cmath>

template<typename T>
struct Node{
	int node_id;
	T x,y,z;
	std::vector<int> link_elem_ids;
	Node(){} //used for trivial initialization: std::vector<Node<T>> nodes_in_elem(elem_type);
	Node(int nid,T x_, T y_, T z_): node_id(nid),x(x_),y(y_),z(z_),link_elem_ids(0) {}
	void print_info(std::ofstream & logfile);
};

struct Edge{
	int edge_id, node1, node2;
};

template<typename T>
class Element {
public:
	int elem_id;
	int attrib;
	std::vector<int> node_ids;
	std::vector<Edge> edges;
	T area;
	Element(int id, int object_id, int node_num): elem_id(id), attrib(object_id), node_ids(node_num), edges(node_num), area(0){}
	virtual T calArea(std::vector<Node<T>> * nodes_ptr) = 0;
	virtual ~Element(){} //The destructor is virtual to ensure that the derived class's destructor is properly called when deleting an object through a base class pointer.
};

template <typename T>
class Triangle : public Element<T> {
public:
	Triangle(int id, int object_id): Element<T>(id,object_id,3){}
	T calArea(std::vector<Node<T>> * nodes_ptr) override;
	~Triangle(){}
};


template<typename T>
class Mesh {
public:
	int num_elem;
	int num_node;
	int elem_type;
	std::vector<Element<T>*> elem_ptrs;
	std::vector<Node<T>> nodes;

	Mesh(int Ne, int Nn, int Et);
	~Mesh(){
		for (int i=0; i<num_elem; ++i){
			delete elem_ptrs[i];
		}
	} //deconstructor for the dynamic memory allocation Element<T>*
	
};

template<typename T>
Mesh<T> * LoadMesh(const std::string meshfile);



#endif
