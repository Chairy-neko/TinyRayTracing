#pragma once
#include <glm/glm.hpp>
#include <string>
#include <iostream>
#include <Eigen/Dense>

using namespace glm;

class Triangle
{
public:
    double calAera();            // calulate area of the triangle
    vec3 findBaryCor(vec3 hitp); // fing the barycenter coordinate of the hitpoint

    Triangle() {}

    vec3 v[3];                         // vertex
    vec3 vn[3];                        // vertex normal
    vec2 vt[3];                        // vertex texture
    vec3 normal = vec3(0.0, 0.0, 0.0); // for caculating hitpoints and distances
    vec3 center = vec3(0.0, 0.0, 0.0); // for sorting
    double area = 0.0;                 // for light sampling
    std::string mtl_name = "";
    bool is_emissive = false; // choose Emissive triangle when they are overlapping
    // int id = 0; // for debug
};