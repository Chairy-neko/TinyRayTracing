#pragma once

#include "ray.h"

#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>
#include <Eigen/Dense>
#include <tinyxml2.h>
#include <opencv2/opencv.hpp>
#include <sstream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <iostream>

using namespace tinyxml2;
using namespace glm;
using namespace std;

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
    string mtl_name = "";
    bool isEmissive = false;
    int id = 0; //for debug
};

class Camera
{
public:
    Camera() {}
    ~Camera() {}

    void SetCamera();
    void Print();
    Ray get_ray(float s, float t);//get ray at (s,t)

    double fovy = 90;
    vec3 eye = vec3(278.0, 273.0, -800.0);
    vec3 lookat = vec3(278.0, 273.0, -799.0);
    vec3 up = vec3(0.0, 1.0, 0.0);
    double aspect_ratio = 1.0;

    vec3 lower_left_corner = vec3(0.0, 0.0, 0.0);
    vec3 horizontal = vec3(0.0, 0.0, 0.0);
    vec3 vertical = vec3(0.0, 0.0, 0.0);
};

class Material
{
public:
    Material() {}
    ~Material() {}
    void readinMap();

    //string name = "";              // name of the material
    vec3 Kd = vec3(0.0, 0.0, 0.0); // diffuse
    vec3 Ks = vec3(0.0, 0.0, 0.0); // specular
    vec3 Tr = vec3(0.0, 0.0, 0.0); // transmittance
    float Ns = 1;                  // shiness, the exponent of phong lobe
    float Ni = 1;                  // the Index of Refraction (IOR) of transparent object
    string map_Kd = "";            // map_Kd

    bool isEmissive = false;
    vec3 radiance = vec3(0, 0, 0);
    double area = 0.0;

    vector<Triangle> triangles;
    cv::Mat img;
    int map_height, map_width;
};

class Light
{
public:
    Light() {}
    Light(string m, vec3 r) : mtl_name(m), radiance(r) {}
    ~Light() {}

    string mtl_name = "";
    vec3 radiance = vec3(0.0, 0.0, 0.0);
};

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
};
