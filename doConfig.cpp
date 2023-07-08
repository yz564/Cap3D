#include "doConfig.h"

template <typename T>
Config<T>::Config(T eps_r, int Nm, std::vector<T> phi, int tds, const std::string mf):\
relative_permitivity(eps_r), metal_num(Nm), electrostatic_potential(phi), allowedthreads(tds), mesh_file(mf){
	Config<T>::print();
}

template <typename T>
void Config<T>::print() {
	std::cout << "[Config] constructor: \n";
	std::cout << "Homogeneous background permitivity : " << relative_permitivity << "F/m \n";
	std::cout << "The number of metal objective: " << metal_num << "\n";
	for (int i=0; i<metal_num; ++i){
		std::cout << "The electrostatic potential of metal objective " << i+1 << ": " << electrostatic_potential[i] << "V \n";
	}
}


template <typename T>
Config<T> * ReadConfig(const std::string inputfile) {
	std::ifstream readfile(inputfile);
	std::string str;
	std::vector<T> tmp(1);
	get_line_strip_comments(readfile, str); //ignore the comments start by '#'
	ReadLine_helper<T>(str, 1, tmp);
	T eps_r = tmp[0];
	
	std::vector<int> tmp_i(1);
	get_line_strip_comments(readfile, str);
	ReadLine_helper<int>(str, 1, tmp_i);
	int Nm = tmp_i[0];
	
	tmp.resize(Nm);
	get_line_strip_comments(readfile, str);
	ReadLine_helper<T>(str, Nm, tmp);
	std::vector<T> phi = tmp;
	
	get_line_strip_comments(readfile, str);
	ReadLine_helper<int>(str, 1, tmp_i);
	int tds = tmp_i[0];
	
	get_line_strip_comments(readfile, str);
	std::string mf=str;
	
	Config<T> * config_ptr = new Config<T>(eps_r, Nm, phi, tds, mf);
	return config_ptr;
}
//explicit instantiations for commonly used types
//otherwise it results in "undefined reference" errors in linker
template Config<float> * ReadConfig<float>(const std::string inputfile);
template Config<double> * ReadConfig<double>(const std::string inputfile);

std::istream& get_line_strip_comments(std::istream& stm, std::string& str)
{
	if (std::getline(stm, str))
	{
		auto pos = str.find("#");
		if (pos == 0) return get_line_strip_comments(stm, str);
		else if (pos != std::string::npos) str.erase(pos);
	}
	return stm;
}

template<typename T>
void ReadLine_helper(const std::string & str, const int num, std::vector<T> & data) {
	std::string word = "";
	int n = 0;
	for (size_t i = 0; i < str.length(); ++i) {
		char x = str[i];
		if (x == '\0' || n >= num) {
			break;
		}
		else if (x != ' ') {
			word += x;
		}
		else if (word.length() > 0) {
			if (std::is_floating_point<T>::value) {
				data[n] = std::stod(word);
			}
			else {
				data[n] = std::stoi(word);
			}
			word = "";
			n++;
		}
	}
	if (n < num && word.length()>0) {
		if (std::is_floating_point<T>::value) {
			data[n] = std::stod(word);
		}
		else {
			data[n] = std::stoi(word);
		}
	}
}
//explicit instantiations for commonly used types
//otherwise it results in "undefined reference" errors in linker
template void ReadLine_helper<int>(const std::string & str, const int num, std::vector<int> & data);
template void ReadLine_helper<float>(const std::string & str, const int num, std::vector<float> & data);
template void ReadLine_helper<double>(const std::string & str, const int num, std::vector<double> & data);

