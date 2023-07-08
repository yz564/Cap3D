#include "doConfig.h"
#include "doMesh.h"
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
	
	//const std::string configfile="config.txt";
	Config<T> * config=ReadConfig<T>("config.txt");
	
	Mesh<T> * mesh = LoadMesh<T>(config->mesh_file);
	
	print_time_cost("Load the mesh takes ", logfile, timer);
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
