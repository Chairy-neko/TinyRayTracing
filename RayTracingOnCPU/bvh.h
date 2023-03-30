#pragma once
#include "scene.h"
#include <algorithm>

const float INF = 114514.0f;

struct HitRecord
{
    bool is_hit = false;
    float distance = INF;
    vec3 hitpoint = vec3(0.0, 0.0, 0.0);
    vec3 direction = vec3(0.0, 0.0, 0.0);
    vec3 pn = vec3(0);
    Triangle triangle;
};
struct BVHNode
{
    BVHNode *left = NULL;
    BVHNode *right = NULL;
    int index = 0, num = 0;
    vec3 AA = vec3(0), BB = vec3(0); // bounding box
};

bool cmpx(const Triangle &t1, const Triangle &t2); // sort triangles according to their center
bool cmpy(const Triangle &t1, const Triangle &t2);
bool cmpz(const Triangle &t1, const Triangle &t2);

BVHNode *buildBVH(std::vector<Triangle> &triangles, int l, int r, int leaf_num);    // build BVH using SAH
HitRecord interactTriangle(Triangle triangle, Ray ray);                             // ray interact with triangle
HitRecord interactBVHNode(Ray ray, std::vector<Triangle> &triangles, int l, int r); // find intercted triangle in BVHNode
float interactAABB(Ray r, vec3 AA, vec3 BB);                                        // interact with AABB, return false if not interact
HitRecord traverseBVH(Ray ray, std::vector<Triangle> &triangles, BVHNode *root);    // traverse BVH