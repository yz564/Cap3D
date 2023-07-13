#include "doConfig.h"
#include "doMesh.h"
#include "doMatrix.h"
#include <ctime>
#include <omp.h>


void print_time_cost(const std::string& message, std::ofstream & logfile, std::vector<clock_t>& timer) {
	std::cout << message << (clock() - timer.back()) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
	logfile << message << (clock() - timer.back()) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
	timer.push_back(clock());
}

template <typename T>
void Cal3dCaps(){
	std::ofstream logfile;
	logfile.open("cap3d.log", std::ios_base::trunc);
	std::vector<clock_t> timer;
	std::cout<<"==================start=================="<<std::endl;
	logfile << "==================start==================" << std::endl;
	timer.push_back(clock());
	
	//configutation
	Config<T> * config=ReadConfig<T>("config.txt");
	omp_set_num_threads(config->allowedthreads);
	Eigen::initParallel();
	//Eigen::setNbThreads(config->allowedthreads);
	//mesh
	Mesh<T> * mesh = LoadMesh<T>(config->mesh_file);
	//basis function
	/* test
	for (int i=0; i<mesh->num_node;++i){
		mesh->nodes[i].print_info(logfile);
	}
	*/
	print_time_cost("Load the mesh takes ", logfile, timer);
	
	BasisFunc<T> * basis = new Basis_0<T>();
	basis->setup(mesh);
	
	Eigen::Matrix<std::complex<T>, -1, -1> * A = CalCoeffMat<T>(basis);
	print_time_cost("Calculate the coefficient matrix takes ", logfile, timer);
	Eigen::Matrix<std::complex<T>, -1, 1> * b = CalRhsVec<T>(basis, &(config->electrostatic_potential));
	print_time_cost("Calculate the rhs vector takes ", logfile, timer);
	
	//Eigen::Matrix<std::complex<T>, -1, 1> x(basis->num_base,1);
	//x.setZero();
	//Eigen::HouseholderQR<Eigen::Matrix<std::complex<T>, -1, -1>> qr(*A);
	//x = qr.solve(*b);
	
	
	
	Eigen::Matrix<std::complex<T>, -1, 1> x;
	
	Eigen::FullPivLU<Eigen::Matrix<std::complex<T>, -1, -1>> lu;
    	lu.setThreshold(0.001); // Set the threshold for partial pivoting
    	lu.compute(*A);
    	int i;
    	#pragma omp parallel for private(i)
    	for (i = 0; i < (*A).rows(); ++i) {
    		Eigen::Matrix<std::complex<T>, -1, 1> bi = (*b).segment(i,1); // Extract the current row of the right-hand side vector
        	x.segment(i,1) = lu.solve(bi);
    	}
	
	//x= (*A).lu().solve(*b);
	//std::cout<<"A.min = "<<(*A).minCoeff()<<std::endl;
	//std::cout<<"A.max = "<<(*A).maxCoeff()<<std::endl;
	std::cout<<"A.sum = "<<(*A).sum()<<std::endl;
	std::cout<<"b.sum = "<<(*b).sum()<<std::endl;
	std::cout<<"x.sum = "<<x.sum()<<std::endl;
	print_time_cost("Solve the matrix equation takes ", logfile, timer);
	std::vector<T> x_(basis->num_base,T(0));
	for(int i=0; i<basis->num_base;++i){
		x_[i]=x(i).real();
	}
	std::vector<T> Q = basis->calQ(x_, config->metal_num);
	for(int i=0; i<config->metal_num;++i){
		std::cout << "Total electric charge on the metal objective " << i+1 << ": " << Q[i] << "C \n";
	}
	
	delete A;
	delete b;
	delete basis;
	delete mesh;
	delete config;
	std::cout << "Total time: " << (clock() - timer.front()) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
	logfile << "Total time: " << (clock() - timer.front()) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
	logfile.close();
}

int main() {
	Cal3dCaps<float>();
	return 0;
}
