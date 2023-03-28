#ifndef MATERIAL_H
#define MATERIAL_H

#include "utils.h"

typedef struct Material
{
    string name = "";               // name of the material
    vec3 Kd = vec3(0.0, 0.0, 0.0);  // diffuse
    vec3 Ks = vec3(0.0, 0.0, 0.0);  // specular
    vec3 Tr = vec3(0.0, 0.0, 0.0);  // transmittance
    float Ns = 1;                   // shiness, the exponent of phong lobe
    float Ni = 1;                   // the Index of Refraction (IOR) of transparent object 
    string map_Kd = "";             // map_Kd

    bool isEmissive = false;      
    vec3 radiance = vec3(0, 0, 0);    
}Material;

typedef struct Light
{
    string mtl_name = "";
    vec3 radiance = vec3(0.0, 0.0, 0.0);
}Light;

#endif