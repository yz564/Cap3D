#ifndef DOMESH_H
#define DOMESH_H

#include <vector>
#include <string>
#include <fstream>
#include <cmath>

template<typename T>
struct Quadrature_point{
	T x,y,z,weight;
	/*
	Quadrature_point& operator=(const Quadrature_point& other){
		if (this == &other){
			return *this;
		}
		x=other.x;
		y=other.y;
		z=other.z;
		weight=other.weight;
		return *this;
	}
	*/
};

struct Base{ //cannot be instantiated
	virtual ~Base(){} //pure virtual function
};

template<typename T>
struct Node: Base{
	int node_id;
	T x,y,z;
	std::vector<int> link_elem_ids;
	Node(){} //used for trivial initialization: std::vector<Node<T>> nodes_in_elem(elem_type);
	Node(int nid,T x, T y, T z): node_id(nid),x(x),y(y),z(z),link_elem_ids(0) {}
	~Node(){}
	void print_info();
};

template<typename T>
struct Edge: Base{
	int edge_id;
	Node<T> node1;
	Node<T> node2;
	~Edge(){}
};

template<typename T>
class Element: public Base{
public:
	int elem_id;
	int attrib;
	std::vector<int> node_ids;
	std::vector<Node<T>*> node_ptrs;
	std::vector<Edge<T>*> edge_ptrs;
	T area;
	Quadrature_point<T> center;
	std::vector<Quadrature_point<T>> qpoints;
	Element(int id, int object_id, int node_num): elem_id(id), attrib(object_id), node_ids(node_num), node_ptrs(node_num), edge_ptrs(node_num), area(0){}
	virtual T calArea() = 0;
	virtual void setQuadrature(int degree) =0;
	virtual ~Element(){} //The destructor is virtual to ensure that the derived class's destructor is properly called when deleting an object through a base class pointer.
};

template <typename T>
class Triangle : public Element<T> {
public:
	Triangle(int id, int object_id): Element<T>(id,object_id,3){}
	T calArea() override;
	void setQuadrature(int degree) override;
	~Triangle(){}
};

template <typename T>
class Quadrilateral : public Element<T> {
public:
	Quadrilateral(int id, int object_id): Element<T>(id,object_id,4){}
	T calArea() override;
	void setQuadrature(int degree) override;
	~Quadrilateral(){}
};


template<typename T>
class Mesh {
public:
	int num_elem;
	int num_node;
	int elem_type;
	int num_attrib;
	std::vector<Element<T>*> elem_ptrs;
	std::vector<Node<T>> nodes;

	Mesh(int Ne, int Nn, int Et, int Na);
	~Mesh(){
		for (int i=0; i<num_elem; ++i){
			delete elem_ptrs[i];
		}
	} //deconstructor for the dynamic memory allocation Element<T>*
	
};

template<typename T>
Mesh<T> * LoadMesh(const std::string meshfile);



#endif
