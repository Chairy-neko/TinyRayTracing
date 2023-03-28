#pragma once
#include "material.h"

//// 光线
//typedef struct Ray
//{
//    vec3 startPoint = vec3(0, 0, 0);    // 起点
//    vec3 direction = vec3(0, 0, 0);     // 方向
//}Ray;

// 光线求交结果
// typedef struct HitResult
// {
//     bool isHit = false;             // 是否命中
//     double distance = 0.0f;         // 与交点的距离
//     vec3 hitPoint = vec3(0, 0, 0);  // 光线命中点
//     Material* material;              // 命中点的表面材质
// }HitResult;

// // 求交结果
// struct HitResult {
//     Triangle* triangle = NULL;
//     float distance = INF;
// };