#include "doMatrix.h"

template <typename T>
T calDistant (Node<T> * p1, Node<T> * p2){
	T dx=p2->x -p1->x;
	T dy=p2->y -p1->y;
	T dz=p2->z -p1->z;
	return sqrt(dx*dx+dy*dy+dz*dz);
}
template float calDistant(Node<float> * p1, Node<float> * p2);
template double calDistant(Node<double> * p1, Node<double> * p2);

template <typename T>
void Basis_0<T>::setup(Mesh<T>* mesh){
	this->mesh=mesh;
	this->num_base=mesh->num_elem;
	this->base.resize(this->num_base);
	for (size_t i=0; i < this->base.size(); ++i){
		this->base[i] = mesh->elem_ptrs[i]->center;
	}
}
template void Basis_0<float>::setup(Mesh<float>* mesh);
template void Basis_0<double>::setup(Mesh<double>* mesh);

template <typename T>
T Basis_0<T>::calDiag(int i){
	T area=this->mesh->elem_ptrs[i]->area;
	return T (1.0/2.0/PERMITIVITY<T>)*sqrt(area/M_PI);
}
template float Basis_0<float>::calDiag(int i);
template double Basis_0<double>::calDiag(int i);


template <typename T>
T Basis_0<T>::calOffDiag(int i, int j){
	T area=this->mesh->elem_ptrs[j]->area;
	Node<T> * r_i = this->mesh->elem_ptrs[i]->center;
	Node<T> * r_j = this->mesh->elem_ptrs[j]->center;
	T dist = calDistant<T>(r_i,r_j);
	return T (1.0/4.0/M_PI/PERMITIVITY<T>)*(area/dist);
}
template float Basis_0<float>::calOffDiag(int i, int j);
template double Basis_0<double>::calOffDiag(int i, int j);


template <typename T>
T Basis_0<T>::calRhs(int i,std::vector<T> * e_potential){
	//if (abs(this->mesh->elem_ptrs[i]->center->z)<=0.001){this->mesh->elem_ptrs[i]->attrib=2;} //temporary fix for the problem mesh
	int attrib=this->mesh->elem_ptrs[i]->attrib;
	return (*e_potential)[attrib-1];
}
template float Basis_0<float>::calRhs(int i,std::vector<float> * e_potential);
template double Basis_0<double>::calRhs(int i,std::vector<double> * e_potential);


template <typename T>
std::vector<T> Basis_0<T>::calQ(std::vector<T>& solution, int num_metals){
	std::vector<T> Q(num_metals,T(0));
	for (int i=0;i< this->num_base; ++i){
		int attrib=this->mesh->elem_ptrs[i]->attrib;
		T area=this->mesh->elem_ptrs[i]->area;
		Q[attrib-1]+=area*solution[i];
	}
	return Q;
}
template std::vector<float> Basis_0<float>::calQ(std::vector<float>& solution, int num_metals);
template std::vector<double> Basis_0<double>::calQ(std::vector<double>& solution, int num_metals);

template <typename T>
Eigen::Matrix<std::complex<T>, -1, -1> * CalCoeffMat(BasisFunc<T> * basis){
	int num_base=basis->num_base;
	Eigen::Matrix<std::complex<T>, -1, -1> * A = new Eigen::Matrix<std::complex<T>, -1, -1>(num_base, num_base);
	//#pragma omp parallel for
	for (int i = 0; i < num_base; ++i) {
		//std::cout << "i="<< i <<" executed by thread " << omp_get_thread_num() << std::endl;
		#pragma omp parallel for
		for (int j = 0; j < num_base; ++j) {
			if (i==j) {
				(*A)(i,j)=basis->calDiag(i);
			}
			else {
				(*A)(i,j)=basis->calOffDiag(i,j);
			}
		}
	}
	return A;
}
template Eigen::Matrix<std::complex<float>, -1, -1> * CalCoeffMat(BasisFunc<float> * basis);
template Eigen::Matrix<std::complex<double>, -1, -1> * CalCoeffMat(BasisFunc<double> * basis);

template <typename T>
Eigen::Matrix<std::complex<T>, -1, -1> * CalRhsVec(BasisFunc<T> * basis, std::vector<T> * e_potential){
	int num_base=basis->num_base;
	Eigen::Matrix<std::complex<T>, -1, -1> * b = new Eigen::Matrix<std::complex<T>, -1, -1>(num_base,1);
	for (int i=0; i<num_base; ++i){
		(*b)(i)=basis->calRhs(i,e_potential);
	}
	return b;
}
template Eigen::Matrix<std::complex<float>, -1, -1> * CalRhsVec(BasisFunc<float> * basis, std::vector<float> * e_potential);
template Eigen::Matrix<std::complex<double>, -1, -1> * CalRhsVec(BasisFunc<double> * basis, std::vector<double> * e_potential);


template <typename T>
Eigen::Matrix<std::complex<T>, -1, -1> * CalRhs(BasisFunc<T> * basis, int num_metal){
	int num_base=basis->num_base;
	Eigen::Matrix<std::complex<T>, -1, -1> * b = new Eigen::Matrix<std::complex<T>, -1, -1>(num_base,num_metal);
	for (int j=0; j<num_metal; ++j){
		std::vector<T> e_potential(num_metal,0);
		e_potential[j]=1;
		for (int i=0; i<num_base; ++i){
			(*b)(i,j)=basis->calRhs(i,&e_potential);
		}
	}
	return b;
}
template Eigen::Matrix<std::complex<float>, -1, -1> * CalRhs(BasisFunc<float> * basis, int num_metal);
template Eigen::Matrix<std::complex<double>, -1, -1> * CalRhs(BasisFunc<double> * basis, int num_metal);
