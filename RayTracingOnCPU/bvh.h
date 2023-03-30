#pragma once
#include "scene.h"
#include <algorithm>

#define INF 114514.0

// 求交结果
struct HitResult
{
    bool isHit = false;
    //bool isInside = false;
    float distance = INF;
    vec3 hitpoint = vec3(0.0, 0.0, 0.0);
    //vec3 normal = vec3(0.0, 0.0, 0.0);
    vec3 viewDir = vec3(0.0, 0.0, 0.0); //击中该点的光线的方向
    vec3 pn = vec3(0);
    //string mtl_name = "";
    //string model_name = "";
    Triangle triangle;
    //int triangle_id = 0;
};

// BVH 树节点
struct BVHNode
{
    BVHNode *left = NULL; // 左右子树索引
    BVHNode *right = NULL;
    int n = 0, index = 0; // 叶子节点信息
    vec3 AA = vec3(0), BB = vec3(0);  // 碰撞盒
};

// 按照三角形中心排序 -- 比较函数
bool cmpx(const Triangle& t1, const Triangle& t2);
bool cmpy(const Triangle& t1, const Triangle& t2);
bool cmpz(const Triangle& t1, const Triangle& t2);

BVHNode* buildBVHwithSAH(std::vector<Triangle>& triangles, int l, int r, int n);// SAH 优化构建 BVH
HitResult hitTriangle(Triangle triangle, Ray ray);// 光线和三角形求交
HitResult hitTriangleArray(Ray ray, std::vector<Triangle>& triangles, int l, int r);// 暴力查数组
float hitAABB(Ray r, vec3 AA, vec3 BB);// 和 aabb 盒子求交，没有交点则返回 -1
HitResult hitBVH(Ray ray, std::vector<Triangle>& triangles, BVHNode* root);// 在 BVH 上遍历求交