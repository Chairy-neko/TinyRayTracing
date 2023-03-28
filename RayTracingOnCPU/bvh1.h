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
bool cmpx(const Triangle &t1, const Triangle &t2)
{
    return t1.center.x < t2.center.x;
}
bool cmpy(const Triangle &t1, const Triangle &t2)
{
    return t1.center.y < t2.center.y;
}
bool cmpz(const Triangle &t1, const Triangle &t2)
{
    return t1.center.z < t2.center.z;
}

// SAH 优化构建 BVH
BVHNode *buildBVHwithSAH(std::vector<Triangle> &triangles, int l, int r, int n)
{
    if (l > r)
        return 0;

    BVHNode *node = new BVHNode();
    node->AA = vec3(1145141919, 1145141919, 1145141919);
    node->BB = vec3(-1145141919, -1145141919, -1145141919);

    // 计算 AABB
    for (int i = l; i <= r; i++)
    {
        // 最小点 AA
        float minx = glm::min(triangles[i].p[0].x, glm::min(triangles[i].p[1].x, triangles[i].p[2].x));
        float miny = glm::min(triangles[i].p[0].y, glm::min(triangles[i].p[1].y, triangles[i].p[2].y));
        float minz = glm::min(triangles[i].p[0].z, glm::min(triangles[i].p[1].z, triangles[i].p[2].z));
        node->AA.x = glm::min(node->AA.x, minx - 0.001f);
        node->AA.y = glm::min(node->AA.y, miny - 0.001f);
        node->AA.z = glm::min(node->AA.z, minz - 0.001f);
        // 最大点 BB
        float maxx = glm::max(triangles[i].p[0].x, glm::max(triangles[i].p[1].x, triangles[i].p[2].x));
        float maxy = glm::max(triangles[i].p[0].y, glm::max(triangles[i].p[1].y, triangles[i].p[2].y));
        float maxz = glm::max(triangles[i].p[0].z, glm::max(triangles[i].p[1].z, triangles[i].p[2].z));
        node->BB.x = glm::max(node->BB.x, maxx + 0.001f);
        node->BB.y = glm::max(node->BB.y, maxy + 0.001f);
        node->BB.z = glm::max(node->BB.z, maxz + 0.001f);
    }

    // 不多于 n 个三角形 返回叶子节点
    if ((r - l + 1) <= n)  {
        node->n = r - l + 1;
        node->index = l;
        return node;
    }

    // 否则递归建树
    float Cost = INF;
    int Axis = 0;
    int Split = (l + r) / 2;
    for (int axis = 0; axis < 3; axis++)
    {
        // 分别按 x，y，z 轴排序
        if (axis == 0)
            std::sort(triangles.begin() + l, triangles.begin() + r + 1, cmpx);
        if (axis == 1)
            std::sort(triangles.begin() + l, triangles.begin() + r + 1, cmpy);
        if (axis == 2)
            std::sort(triangles.begin() + l, triangles.begin() + r + 1, cmpz);

        // leftMax[i]: [l, i] 中最大的 xyz 值
        // leftMin[i]: [l, i] 中最小的 xyz 值
        std::vector<vec3> leftMax(r - l + 1, vec3(-INF, -INF, -INF));
        std::vector<vec3> leftMin(r - l + 1, vec3(INF, INF, INF));
        // 计算前缀 注意 i-l 以对齐到下标 0
        for (int i = l; i <= r; i++)
        {
            Triangle t = triangles[i];
            int bias = (i == l) ? 0 : 1; // 第一个元素特殊处理

            leftMax[i - l].x = glm::max(leftMax[i - l - bias].x, glm::max(t.p[0].x, glm::max(t.p[1].x, t.p[2].x)));
            leftMax[i - l].y = glm::max(leftMax[i - l - bias].y, glm::max(t.p[0].y, glm::max(t.p[1].y, t.p[2].y)));
            leftMax[i - l].z = glm::max(leftMax[i - l - bias].z, glm::max(t.p[0].z, glm::max(t.p[1].z, t.p[2].z)));

            leftMin[i - l].x = glm::min(leftMin[i - l - bias].x, glm::min(t.p[0].x, glm::min(t.p[1].x, t.p[2].x)));
            leftMin[i - l].y = glm::min(leftMin[i - l - bias].y, glm::min(t.p[0].y, glm::min(t.p[1].y, t.p[2].y)));
            leftMin[i - l].z = glm::min(leftMin[i - l - bias].z, glm::min(t.p[0].z, glm::min(t.p[1].z, t.p[2].z)));
        }

        // rightMax[i]: [i, r] 中最大的 xyz 值
        // rightMin[i]: [i, r] 中最小的 xyz 值
        std::vector<vec3> rightMax(r - l + 1, vec3(-INF, -INF, -INF));
        std::vector<vec3> rightMin(r - l + 1, vec3(INF, INF, INF));
        // 计算后缀 注意 i-l 以对齐到下标 0
        for (int i = r; i >= l; i--)
        {
            Triangle t = triangles[i];
            int bias = (i == r) ? 0 : 1; // 第一个元素特殊处理

            rightMax[i - l].x = glm::max(rightMax[i - l + bias].x, glm::max(t.p[0].x, glm::max(t.p[1].x, t.p[2].x)));
            rightMax[i - l].y = glm::max(rightMax[i - l + bias].y, glm::max(t.p[0].y, glm::max(t.p[1].y, t.p[2].y)));
            rightMax[i - l].z = glm::max(rightMax[i - l + bias].z, glm::max(t.p[0].z, glm::max(t.p[1].z, t.p[2].z)));

            rightMin[i - l].x = glm::min(rightMin[i - l + bias].x, glm::min(t.p[0].x, glm::min(t.p[1].x, t.p[2].x)));
            rightMin[i - l].y = glm::min(rightMin[i - l + bias].y, glm::min(t.p[0].y, glm::min(t.p[1].y, t.p[2].y)));
            rightMin[i - l].z = glm::min(rightMin[i - l + bias].z, glm::min(t.p[0].z, glm::min(t.p[1].z, t.p[2].z)));
        }

        // 遍历寻找分割
        float cost = INF;
        int split = l;
        for (int i = l; i <= r - 1; i++)
        {
            float lenx, leny, lenz;
            // 左侧 [l, i]
            vec3 leftAA = leftMin[i - l];
            vec3 leftBB = leftMax[i - l];
            lenx = leftBB.x - leftAA.x;
            leny = leftBB.y - leftAA.y;
            lenz = leftBB.z - leftAA.z;
            float leftS = 2.0 * ((lenx * leny) + (lenx * lenz) + (leny * lenz));
            float leftCost = leftS * (i - l + 1);

            // 右侧 [i+1, r]
            vec3 rightAA = rightMin[i + 1 - l];
            vec3 rightBB = rightMax[i + 1 - l];
            lenx = rightBB.x - rightAA.x;
            leny = rightBB.y - rightAA.y;
            lenz = rightBB.z - rightAA.z;
            float rightS = 2.0 * ((lenx * leny) + (lenx * lenz) + (leny * lenz));
            float rightCost = rightS * (r - i);

            // 记录每个分割的最小答案
            float totalCost = leftCost + rightCost;
            if (totalCost < cost)
            {
                cost = totalCost;
                split = i;
            }
        }
        // 记录每个轴的最佳答案
        if (cost < Cost)
        {
            Cost = cost;
            Axis = axis;
            Split = split;
        }
    }

    // 按最佳轴分割
    if (Axis == 0)
        std::sort(triangles.begin() + l, triangles.begin() + r + 1, cmpx);
    if (Axis == 1)
        std::sort(triangles.begin() + l, triangles.begin() + r + 1, cmpy);
    if (Axis == 2)
        std::sort(triangles.begin() + l, triangles.begin() + r + 1, cmpz);

    // 递归
    node->left = buildBVHwithSAH(triangles, l, Split, n);
    node->right = buildBVHwithSAH(triangles, Split + 1, r, n);

    return node;
}

// 光线和三角形求交
HitResult hitTriangle(Triangle triangle, Ray ray)
{
    HitResult res;
    vec3 p1 = triangle.p[0], p2 = triangle.p[1], p3 = triangle.p[2];
    vec3 S = ray.startPoint;                     // 射线起点
    vec3 d = ray.direction;                      // 射线方向
    vec3 N = triangle.normal;
    //vec3 N = normalize(cross(p2 - p1, p3 - p1)); // 法向量
    // if (dot(N, d) > 0.0f)
    // {
    //     N = -N; // 获取正确的法向量
    //     //res.isInside = true;
    // }
    // 如果视线和三角形平行
    if (fabs(dot(N, d)) < 0.00001f)
        return res;

    // 距离
    float t = (dot(p1-S, N)) / dot(d, N);
    if (t < 0.0005f)
        return res; // 如果三角形在光线背面

    // 交点计算
    vec3 P = S + d * t;

    // 判断交点是否在三角形中
    vec3 c1 = cross(p2 - p1, P - p1);
    vec3 c2 = cross(p3 - p2, P - p2);
    vec3 c3 = cross(p1 - p3, P - p3);
    double dir1 = dot(c1, N), dir2 = dot(c2, N), dir3 = dot(c3, N);
    bool r1 = dir1 > 0 && dir2 > 0 && dir3 > 0;
    bool r2 = dir1 < 0 && dir2 < 0 && dir3 < 0;
    /*bool r1 = (dot(c1, N) > 0 && dot(c2, N) > 0 && dot(c3, N) > 0);
    bool r2 = (dot(c1, N) < 0 && dot(c2, N) < 0 && dot(c3, N) < 0);*/

    // 命中，封装返回结果
    if (r1 || r2) {
        // if (triangle.mtl_name == "back:LeftWall")
        //     cout << "triangle id: " << triangle.id << endl;
        res.isHit = true;
        res.hitpoint = P;
        res.distance = t;
        res.viewDir = d;
    }

    return res;
}

// 暴力查数组
HitResult hitTriangleArray(Ray ray, std::vector<Triangle> &triangles, int l, int r)
{
    HitResult res;
    for (int i = l; i <= r; i++)
    {
        HitResult restmp = hitTriangle(triangles[i], ray);
        if(restmp.isHit && restmp.distance < res.distance){
            res = restmp;
            res.triangle = triangles[i];
            vec3 barycenter = res.triangle.findBaryCor(res.hitpoint);
            res.pn = (res.triangle.vn[0] * barycenter.x) + (res.triangle.vn[1] * barycenter.y) + (res.triangle.vn[2] * barycenter.z);
        }
    }
    return res;
}

// 和 aabb 盒子求交，没有交点则返回 -1
float hitAABB(Ray r, vec3 AA, vec3 BB)
{
    // 1.0 / direction
    vec3 invdir = vec3(1.0 / r.direction.x, 1.0 / r.direction.y, 1.0 / r.direction.z);

    vec3 in = (BB - r.startPoint) * invdir;
    vec3 out = (AA - r.startPoint) * invdir;

    vec3 tmax = glm::max(in, out);
    vec3 tmin = glm::min(in, out);

    float t1 = glm::min(tmax.x, glm::min(tmax.y, tmax.z));
    float t0 = glm::max(tmin.x, glm::max(tmin.y, tmin.z));

    return (t1 >= t0) ? ((t0 > 0.0) ? (t0) : (t1)) : (-1);
}

// 在 BVH 上遍历求交
HitResult hitBVH(Ray ray, std::vector<Triangle> &triangles, BVHNode *root)
{
    if (root == NULL)
        return HitResult();

    // 是叶子 暴力查
    if (root->n > 0)
    {
        return hitTriangleArray(ray, triangles, root->index, root->n + root->index - 1);
    }

    // 和左右子树 AABB 求交
    float d1 = INF, d2 = INF;
    if (root->left)
        d1 = hitAABB(ray, root->left->AA, root->left->BB);
    if (root->right)
        d2 = hitAABB(ray, root->right->AA, root->right->BB);

    // 递归结果
    HitResult r1, r2;
    if (d1 > 0)
        r1 = hitBVH(ray, triangles, root->left);
    if (d2 > 0)
        r2 = hitBVH(ray, triangles, root->right);
    //if r1.d == r2.d, return light
    if(r1.distance == r2.distance)
        if(r1.triangle.isEmissive)
            return r1;
        else
            return r2;

    return r1.distance < r2.distance ? r1 : r2;
}