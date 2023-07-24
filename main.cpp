#include "doConfig.h"
#include "doMesh.h"
#include "doMatrix.h"
#include <ctime>
#include <omp.h>
#include <mkl.h>
#include <cstdlib>
#include <cassert>

#define EIGEN_USE_MKL_ALL

std::ofstream logfile("cap3d.log");

void print_time_cost(const std::string& message, std::vector<clock_t>& timer) {
	std::cout << message << (clock() - timer.back()) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
	logfile << message << (clock() - timer.back()) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
	timer.push_back(clock());
}

template<typename T>
void printTransposeMatrix(const std::vector<std::vector<T>>& matrix)
{
    for (size_t j = 0; j < matrix[0].size(); ++j)
    {
        for (size_t i = 0; i < matrix.size(); ++i)
        {
            std::cout << matrix[i][j] << "  ";
            logfile << matrix[i][j] << "  ";
        }
        std::cout << std::endl;
        logfile << std::endl;
    }
}

template<typename T>
void get_cap_matrix(int metal_num, BasisFunc<T> * basis, Eigen::Matrix<std::complex<T>, -1, -1> * b){
	std::vector<std::vector<T>> cap_mat(metal_num);
	for(int j=0; j<metal_num;++j){
		std::vector<T> x(basis->num_base,T(0));
		for(int i=0; i<basis->num_base;++i){
			x[i]=(*b)(i,j).real();
		}
		std::vector<T> Q = basis->calQ(x, metal_num);
		cap_mat[j]=Q;
	}
	std::cout << "The capacitance matrix:"<<std::endl;
	logfile << "The capacitance matrix:"<<std::endl;
	printTransposeMatrix(cap_mat);
}

template<typename T>
void get_charge_distribution(int metal_num, BasisFunc<T> * basis, Eigen::Matrix<std::complex<T>, -1, -1> * b){
	std::vector<T> x(basis->num_base,T(0));
	for(int i=0; i<basis->num_base;++i){
		x[i]=(*b)(i,0).real();
	}
	std::vector<T> Q = basis->calQ(x, metal_num);
	for(int i=0; i<metal_num;++i){
		std::cout << "Total electric charge on the metal objective " << i+1 << ": " << Q[i] << "C \n";
		logfile << "Total electric charge on the metal objective " << i+1 << ": " << Q[i] << "C \n";
	}
	
	//export results to plot
	std::ofstream outFile("mesh_charge_density.csv");
	outFile << "coordinates" << std::endl;
	outFile << "x,y,z" << std::endl;
    	for (int i=0; i<basis->mesh->num_node; ++i) {
        	outFile << basis->mesh->nodes[i].x << "," << basis->mesh->nodes[i].y << "," << basis->mesh->nodes[i].z << std::endl;
    	}
    	outFile << std::endl; //to separate the tables
    	outFile << "elements" << std::endl;
	for (int j=0; j<basis->mesh->elem_type-1; ++j){
    		outFile << "vertex"<<j+1<<",";	
    	}
    	outFile << "vertex"<<basis->mesh->elem_type<<std::endl;
    	for (int i=0; i<basis->mesh->num_elem; ++i) {
    		for (int j=0; j<basis->mesh->elem_type-1; ++j){
    			outFile << basis->mesh->elem_ptrs[i]->node_ids[j]<<",";
    		}
        	outFile << basis->mesh->elem_ptrs[i]->node_ids[basis->mesh->elem_type-1] <<std::endl;
    	}
    	outFile << std::endl;
    	outFile << "charge_density" << std::endl;
	outFile << "q" << std::endl;
	for (int i=0; i<basis->num_base; ++i) {
        	outFile << x[i] << std::endl;
    	}
    	outFile << std::endl;
    	outFile.close();
}

template <typename T>
void Cal3dCaps(){
	std::vector<clock_t> timer;
	std::cout<<"==================start=================="<<std::endl;
	logfile << "==================start==================" << std::endl;
	timer.push_back(clock());
	
	//configutation
	Config<T> * config=ReadConfig<T>("config.txt");
	const int numThreads = omp_get_max_threads();
	omp_set_num_threads(numThreads);
	mkl_set_num_threads(numThreads);
	std::cout << "Thread numbers: " << numThreads << std::endl;
	Eigen::initParallel();
	//Eigen::setNbThreads(config->allowedthreads);
	//mesh
	Mesh<T> * mesh = LoadMesh<T>(config->mesh_file);
	print_time_cost("Load the mesh takes ", timer);
	
	//basis function
	BasisFunc<T> * basis = new Basis_0<T>(); //0-order basis function and delta testing function
	basis->setup(mesh);
	
	bool is_cap_matrix = std::getenv("CAP_MATRIX");
	Eigen::Matrix<std::complex<T>, -1, -1> * b;
	int metal_num = mesh->num_attrib;
	if (is_cap_matrix){
		b = CalRhs<T>(basis, metal_num); //identity matrix on the right hand side
		
	}else{
		b = CalRhsVec<T>(basis, &(config->electrostatic_potential));
		try{
			assert(metal_num==config->metal_num);
		}catch(const std::exception& e){
			std::cerr << "Error [check metal_num in config.txt]: " << e.what() << std::endl;
		}
	}
	print_time_cost("Calculate the rhs vector takes ", timer);
	
	//prepare the coefficient matrix of double integrations
	Eigen::Matrix<std::complex<T>, -1, -1> * A = CalCoeffMat<T>(basis);
	print_time_cost("Calculate the coefficient matrix takes ", timer);
	
	const int n= basis->num_base;
	//Eigen::Matrix<std::complex<T>, -1, 1> x(n);
	//x= (*A).lu().solve(*b);
	
	// Solve the dense matrix equation Ax = b using MKL
	int nrhs = b->cols();
	int lda = A->rows();
	int ldb = b->rows();
	int info;
	int* ipiv = new int[n];
	

	if(std::is_same<T, float>::value){
		MKL_Complex8* mkl_A = (MKL_Complex8*)mkl_malloc(n * n * sizeof(MKL_Complex8), 64); //Memory Alignment: For better performance
		std::copy(A->data(), A->data() + (n * n), reinterpret_cast<std::complex<float>*>(mkl_A));
		delete A;
		cgesv(&n, &nrhs, mkl_A, &lda, ipiv, reinterpret_cast<MKL_Complex8*>(b->data()), &ldb, &info);
		delete mkl_A;
	}else{
		MKL_Complex16* mkl_A = (MKL_Complex16*)mkl_malloc(n * n * sizeof(MKL_Complex16), 64);
		std::copy(A->data(), A->data() + (n * n), reinterpret_cast<std::complex<double>*>(mkl_A));
		delete A;
		zgesv(&n, &nrhs, mkl_A, &lda, ipiv, reinterpret_cast<MKL_Complex16*>(b->data()), &ldb, &info);
		delete mkl_A;
	}

	if (info > 0) {
		std::cout << "Failed to solve the matrix equation." << std::endl;
	}
	print_time_cost("Solve the matrix equation takes ", timer);
	
	//post processing
	if (is_cap_matrix){
		get_cap_matrix(metal_num, basis, b);
	}else{
		get_charge_distribution(metal_num, basis, b);
	}
	print_time_cost("Post processing takes ", timer);

	delete[] ipiv;
	delete b;
	delete basis;
	delete mesh;
	delete config;
	std::cout << "Total time: " << (clock() - timer.front()) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
	logfile << "Total time: " << (clock() - timer.front()) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
	std::cout<<"==================finish=================="<<std::endl;
	logfile << "==================finish==================" << std::endl;
	logfile.close();
}

int main() {
	Cal3dCaps<float>();
	return 0;
}
