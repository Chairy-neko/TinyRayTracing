#pragma once
#include <glm/glm.hpp>
#include <string>
#include <iostream>
#include <Eigen/Dense>

using namespace glm;

class Triangle
{
public:
    double calAera(); // calulate area of the triangle
    vec3 findBaryCor(vec3 hitp); // fing the barycenter coordinate of the triangle

    Triangle() {}
    vec3 p[3];
    vec3 vn[3];
    vec2 vt[3];
    vec3 normal = vec3(0.0, 0.0, 0.0);
    vec3 center = vec3(0.0, 0.0, 0.0);
    double area = 0.0;
    std::string mtl_name = "";
    bool isEmissive = false;
    int id = 0; //for debug
};