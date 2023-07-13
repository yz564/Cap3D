#ifndef DOMATRIX_H
#define DOMATRIX_H

#include "doMesh.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/IterativeSolvers>
#include <math.h>
#include <cmath>
#include <omp.h>

template<typename T>
constexpr T PERMITIVITY = static_cast<T>(8.854e-12);

template <typename T>
struct BasisFunc{
	int num_base;
	Mesh<T>* mesh;
	std::vector<Base *> base;
	virtual void setup(Mesh<T>* mesh)=0;
	virtual T calDiag(int i)=0;
	virtual T calOffDiag(int i, int j)=0;
	virtual T calRhs(int i,std::vector<T> * e_potential)=0;
	virtual std::vector<T> calQ(std::vector<T> & solution,int num_metals)=0;
	virtual ~BasisFunc(){}
};

template <typename T>
struct Basis_0: public BasisFunc<T>{
	void setup(Mesh<T>* mesh) override;
	T calDiag(int i);
	T calOffDiag(int i, int j);
	T calRhs(int i,std::vector<T> * e_potential);
	std::vector<T> calQ(std::vector<T> & solution,int num_metals);
	~Basis_0(){} // should not specify deconstructor because Base * is deleted in Triangle
};

template <typename T>
Eigen::Matrix<std::complex<T>, -1, -1> * CalCoeffMat(BasisFunc<T> * basis);

template <typename T>
Eigen::Matrix<std::complex<T>, -1, 1> * CalRhsVec(BasisFunc<T> * basis, std::vector<T> * e_potential);

template <typename T>
T calDistant (Node<T> * p1, Node<T> * p2);

#endif
