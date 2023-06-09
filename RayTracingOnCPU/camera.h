#pragma once
#include "ray.h"
#include <glm/glm.hpp>

class Camera
{
public:
    Camera() {}
    ~Camera() {}

    void setCamera();// set camera 
    Ray getRay(float u, float v);// get ray at (u,v)
    void Print();// print info for debugging

    double fovy = 90;
    vec3 eye = vec3(278.0, 273.0, -800.0);
    vec3 lookat = vec3(278.0, 273.0, -799.0);
    vec3 up = vec3(0.0, 1.0, 0.0);
    double aspect_ratio = 1.0;

    vec3 lower_left_corner = vec3(0.0, 0.0, 0.0);
    vec3 horizontal = vec3(0.0, 0.0, 0.0);
    vec3 vertical = vec3(0.0, 0.0, 0.0);
};
