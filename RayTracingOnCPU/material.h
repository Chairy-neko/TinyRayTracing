#pragma once
#include "triangle.h"
#include <glm/glm.hpp>
#include <string>
#include <vector>
#include <opencv2/opencv.hpp>

using namespace glm;
using namespace std;

class Material
{
public:
    Material() {}
    ~Material() {}
    void readinMap(); //for texture

    vec3 Kd = vec3(0.0, 0.0, 0.0); // diffuse
    vec3 Ks = vec3(0.0, 0.0, 0.0); // specular
    vec3 Tr = vec3(0.0, 0.0, 0.0); // transmittance
    float Ns = 1;                  // shiness, the exponent of phong lobe
    float Ni = 1;                  // the Index of Refraction (IOR) of transparent object
    string map_Kd = "";            // map_Kd

    bool is_emissive = false;
    vec3 radiance = vec3(0, 0, 0);
    double area = 0.0;

    vector<Triangle> triangles; // for light sampling

    cv::Mat img;// for texture
    int map_height, map_width;
};