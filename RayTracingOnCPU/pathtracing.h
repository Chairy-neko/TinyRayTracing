#pragma once
#include "bvh.h"

#include <random>
#include <glm/glm.hpp>
#include <iostream>

using namespace std;
using namespace glm;

const float PI = 3.1415926f;
const float P_RR = 0.8; // russian roulette probability

bool RR(double pr);                              // Russian Roulette
vec3 Sample(vec3 &direction, int type, double Ns); // importance sampling
Ray nextRay(HitRecord &res, vec3 dir, Scene &scene);
vec3 shade(HitRecord &res, vec3 dir, Scene &scene, BVHNode *root);