#include "doMesh.h"
#include <assert.h>

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
	int num_blk = tmp[0];
	int num_elem = tmp[1];
	int num_node = tmp[2];
	get_line_strip_comments(readfile, str); //third line
	ReadLine_helper<int>(str, 1, tmp);
	int elem_type = tmp[0];
	Mesh<T> * ans = new Mesh<T>(num_elem, num_node, elem_type);
	int count = 0;
	ans->elems.resize(num_elem);
	ans->attribs.resize(num_elem);
	ans->nodes.resize(num_node);
	for (int i = 0; i < num_blk; ++i) {
		get_line_strip_comments(readfile, str);
		ReadLine_helper<int>(str, 3, tmp);
		int num_elem_blk = tmp[0];
		//int blk_id = tmp[1];
		int attrib_blk = tmp[2];
		std::vector<int> data(elem_type);
		for (int j = 0; j < num_elem_blk; ++j) {
			std::getline(readfile, str);
			ReadLine_helper<int>(str, elem_type, data);
			ans->elems[count]=data; //the meshfile is 1-based
			ans->attribs[count]=attrib_blk;
			count++;
		}
	}
	assert(count == num_elem);
	//start to read the node coordinates (three-dimensional)
	std::vector<T> data2(3);
	for (int i = 0; i < num_node; ++i) {
		std::getline(readfile, str);
		ReadLine_helper<T>(str, 3, data2);
		ans->nodes[i]=data2;
	}
	return ans;
}
template Mesh<float> * LoadMesh<float>(const std::string meshfile);
template Mesh<double> * LoadMesh<double>(const std::string meshfile);

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


