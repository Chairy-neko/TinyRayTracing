#pragma once

//#include "triangle.h"
#include "camera.h"
#include "material.h"
#include "light.h"
//#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>
//#include <Eigen/Dense>
#include <tinyxml2.h>
//#include <opencv2/opencv.hpp>
#include <sstream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <iostream>

using namespace tinyxml2;
using namespace glm;
using namespace std;

class Scene
{
public:
    Scene(/* args */) {}
    ~Scene() {}

    void readxml(string xml_path); // get camera & lights
    void readmtl(string mtl_path, string base_dir); // get materials
    void readobj(string obj_path); // get triangles

    int img_width;
    int img_height;
    vector<Triangle> triangles; // obj models
    vector<Light> lights;
    unordered_map<string, Material> materials;
    Camera camera;

    double total_area = 0.0;
    vector<Triangle> light_triangles;
};
