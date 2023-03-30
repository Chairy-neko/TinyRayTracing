#pragma once
#include <glm/glm.hpp>
#include <string>
#include <iostream>

using namespace std;
using namespace glm;

class Light
{
public:
    Light() {}
    Light(string m, vec3 r) : mtl_name(m), radiance(r) {}
    ~Light() {}

    string mtl_name = "";
    vec3 radiance = vec3(0.0, 0.0, 0.0);
};