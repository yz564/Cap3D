#include "doConfig.h"
#include "doMesh.h"
#include "doMatrix.h"
#include <ctime>


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
	//mesh
	Mesh<T> * mesh = LoadMesh<T>(config->mesh_file);
	//basis function
	/* test
	for (int i=0; i<mesh->num_node;++i){
		mesh->nodes[i].print_info(logfile);
	}
	*/
	BasisFunc<T> * basis = new Basis_0<T>();
	basis->setup(mesh);
	
	Eigen::Matrix<std::complex<T>, -1, -1> * A = CalCoeffMat<T>(basis);
	
	Eigen::Matrix<std::complex<T>, -1, -1> * b = CalRhsVec<T>(basis, &(config->electrostatic_potential));
	
	//Eigen::Matrix<std::complex<T>, -1, -1> x(basis->num_base,1);
	//x= (*A).lu().solve(*b);
	//std::cout<<"A.min = "<<(*A).minCoeff()<<std::endl;
	//std::cout<<"A.max = "<<(*A).maxCoeff()<<std::endl;
	std::cout<<"A.sum = "<<(*A).sum()<<std::endl;
	std::cout<<"b.sum = "<<(*b).sum()<<std::endl;
	//std::cout<<"x.sum = "<<x.sum()<<std::endl;
	print_time_cost("Load the mesh takes ", logfile, timer);
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
	Cal3dCaps<double>();
	return 0;
}
