#pragma once
#include <glm/glm.hpp>
using namespace glm;

const int DIFFUSE      = 0;
const int SPECULAR     = 1;
const int TRANSMISSION = 2;
const int INVALID = 3;

class Ray
{
public:
    Ray() {}
    Ray(vec3 s, vec3 d) :startpoint(s), direction(d) {}
    Ray(vec3 s, vec3 d, int r) :startpoint(s), direction(d), ray_type(r) {}
    ~Ray() {}

    vec3 startpoint = vec3(0.0, 0.0, 0.0);
    vec3 direction = vec3(0.0, 0.0, 0.0);
    int ray_type = INVALID;
};