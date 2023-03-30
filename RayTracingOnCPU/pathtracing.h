#pragma once
#include "bvh.h"

#include <random>
#include <glm/glm.hpp>
#include <iostream>

using namespace std;
using namespace glm;

const float PI = 3.1415926f;
const float P_RR = 0.8; // russian roulette probability

extern int num_kd;
extern int num_ks; // for debugging purposes
extern int dark;

bool russian_Roulette(double pr);
vec3 BRDFImportanceSampling(vec3& direction, int type, double Ns);// importance sampling based on BRDF
Ray nextRay(HitResult& res, vec3 dir, Scene& scene);
vec3 shade(HitResult& res, vec3 dir, Scene& scene, BVHNode* root);
vec3 directIllumination(HitResult &res, vec3 dir, Scene &scene);