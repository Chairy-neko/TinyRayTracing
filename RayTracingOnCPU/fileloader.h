#pragma once
#include <tinyxml2.h>
#include <string>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <sstream>
#include "camera.h"
#include "material.h"
#include "model.h"
using namespace tinyxml2;
using namespace std;

void readxml(string filepath,
			 Camera &camera,
			 vector<Light> &lights,
			 int &img_width,
			 int &img_height)
{
	XMLDocument doc;
	XMLError error = doc.LoadFile(filepath.c_str());
	if (error != XML_SUCCESS)
	{
		std::cout << "Read xml failed." << doc.ErrorStr() << std::endl;
		exit(1);
	}
	XMLElement *xml_camera = doc.RootElement();
	img_width = stoi(xml_camera->Attribute("width"));
	img_height = stoi(xml_camera->Attribute("height"));
	camera.fovy = stof(xml_camera->Attribute("fovy"));
	XMLElement *xml_eye = xml_camera->FirstChildElement("eye");
	camera.eye = vec3(stof(xml_eye->Attribute("x")), stof(xml_eye->Attribute("y")), stof(xml_eye->Attribute("z")));
	XMLElement *xml_lookat = xml_camera->FirstChildElement("lookat");
	camera.lookat = vec3(stof(xml_lookat->Attribute("x")), stof(xml_lookat->Attribute("y")), stof(xml_lookat->Attribute("z")));
	XMLElement *xml_up = xml_camera->FirstChildElement("up");
	camera.up = vec3(stof(xml_up->Attribute("x")), stof(xml_up->Attribute("y")), stof(xml_up->Attribute("z")));
	XMLElement *xml_light = xml_camera->NextSiblingElement("light");
	while (xml_light)
	{
		Light light;
		light.mtl_name = xml_light->Attribute("mtlname");
		string radiance_str = xml_light->Attribute("radiance");
		vec3 radiance = vec3(0.0, 0.0, 0.0);
		int tmp = 0, flag = 0;
		for (int i = 0; i < radiance_str.length(); i++)
		{
			if (radiance_str[i] == ',' && flag == 0)
			{
				radiance.x = stof(radiance_str.substr(tmp, i - tmp));
				tmp = i + 1;
				flag = 1;
			}
			else if (radiance_str[i] == ',' && flag == 1)
			{
				radiance.y = stof(radiance_str.substr(tmp, i - tmp));
				tmp = i + 1;
				flag = 2;
			}
			else if (i == radiance_str.length() - 1)
			{
				radiance.z = stof(radiance_str.substr(tmp, i - tmp + 1));
			}
		}
		light.radiance = radiance;
		lights.push_back(light);
		xml_light = xml_light->NextSiblingElement("light");
	}

	// std::cout << "------------------------------------------Read XML Done.------------------------------------------" << std::endl;
}

void readobj(
	/*std::string filepath,
	std::vector<glm::vec3>&points,
	std::vector<glm::vec2>&texcoords,
	std::vector<glm::vec3>&normals*/
	string objfile,
	string mtlfile,
	unordered_map<string, Material> &mmap,
	vector<Triangle> &triangles)
{
	// get mtl file
	// ���ļ���
	std::ifstream fin(mtlfile);
	std::string line;
	if (!fin.is_open())
	{
		std::cout << "Read " << mtlfile << " failed." << std::endl;
		exit(-1);
	}
	string mtlname;
	while (std::getline(fin, line))
	{
		std::istringstream sin(line); // ��һ�е�������Ϊ string stream �������Ҷ�ȡ
		std::string type;
		float x, y, z;
		float n;

		// ��ȡobj�ļ�
		sin >> type;
		if (type == "newmtl")
		{
			Material material;
			sin >> mtlname;
			material.name = mtlname;
			mmap.insert({mtlname, material});
		}
		else if (type == "Kd")
		{
			sin >> x >> y >> z;
			mmap[mtlname].Kd = vec3(x, y, z);
		}
		else if (type == "Ks")
		{
			sin >> x >> y >> z;
			mmap[mtlname].Ks = vec3(x, y, z);
		}
		else if (type == "Tr")
		{
			sin >> x >> y >> z;
			mmap[mtlname].Tr = vec3(x, y, z);
		}
		else if (type == "Ns")
		{
			sin >> n;
			mmap[mtlname].Ns = n;
		}
		else if (type == "Ni")
		{
			sin >> n;
			mmap[mtlname].Ni = n;
		}
	}
	cout << "num of materials  : " << mmap.size() << endl;
	//cout << "------------------------------------------Read MTL Done.------------------------------------------" << endl;

	// get obj file
	// ���ļ���
	std::ifstream fin2(objfile);
	std::string line2;
	if (!fin2.is_open())
	{
		std::cout << "Read " << objfile << " failed." << std::endl;
		exit(-1);
	}

	// ������ȡ
	// int offset = vertices.size();
	// ���ж�ȡ
	std::vector<glm::vec3> vertices;
	std::string mtl_name = "";
	while (std::getline(fin2, line2))
	{
		std::istringstream sin(line2); // ��һ�е�������Ϊ string stream �������Ҷ�ȡ
		std::string type;
		float x, y, z;
		string f[3];

		// ��ȡobj�ļ�
		sin >> type;
		if (type == "v")
		{
			sin >> x >> y >> z;
			vertices.push_back(vec3(x, y, z));
		}
		else if (type == "usemtl")
		{
			sin >> mtl_name;
		}
		else if (type == "f")
		{
			sin >> f[0] >> f[1] >> f[2];
			int tmpindex = 0;
			Triangle triangle;
			for (int fi = 0; fi < 3; fi++)
			{
				for (int index = 0; index < f[fi].length(); index++)
				{
					if (f[fi][index] == '/')
					{
						tmpindex = stoi(f[fi].substr(0, index)) - 1;
						triangle.p[fi] = vertices[tmpindex];
						break;
					}
				}
			}
			triangle.normal = normalize(cross(triangle.p[1] - triangle.p[0], triangle.p[2] - triangle.p[0]));
			triangle.center = (triangle.p[0] + triangle.p[1] + triangle.p[2]) / vec3(3.0, 3.0, 3.0);
			triangle.mtl_name = mtl_name;
			triangles.push_back(triangle);
		}
	}
	cout << "num of vertices  : " << vertices.size() << endl;
	cout << "num of triangles : " << triangles.size() << endl;
	// cout << "------------------------------------------Read OBJ Done.------------------------------------------" << endl;
}
