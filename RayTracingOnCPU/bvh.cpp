#include "bvh.h"

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

BVHNode *buildBVH(std::vector<Triangle> &triangles, int left_index, int right_index, int leaf_num)
{
    if (left_index > right_index)
        return 0;

    BVHNode *node = new BVHNode();
    node->AA = vec3(1145141919, 1145141919, 1145141919);
    node->BB = vec3(-1145141919, -1145141919, -1145141919);
    // caculate AABB
    for (int i = left_index; i <= right_index; i++)
    {
        // minimum point AA
        float x_min = glm::min(triangles[i].v[0].x, glm::min(triangles[i].v[1].x, triangles[i].v[2].x));
        float y_min = glm::min(triangles[i].v[0].y, glm::min(triangles[i].v[1].y, triangles[i].v[2].y));
        float z_min = glm::min(triangles[i].v[0].z, glm::min(triangles[i].v[1].z, triangles[i].v[2].z));
        node->AA.x = glm::min(node->AA.x, x_min - 0.001f);
        node->AA.y = glm::min(node->AA.y, y_min - 0.001f);
        node->AA.z = glm::min(node->AA.z, z_min - 0.001f);
        // maximum point BB
        float x_max = glm::max(triangles[i].v[0].x, glm::max(triangles[i].v[1].x, triangles[i].v[2].x));
        float y_max = glm::max(triangles[i].v[0].y, glm::max(triangles[i].v[1].y, triangles[i].v[2].y));
        float z_max = glm::max(triangles[i].v[0].z, glm::max(triangles[i].v[1].z, triangles[i].v[2].z));
        node->BB.x = glm::max(node->BB.x, x_max + 0.001f);
        node->BB.y = glm::max(node->BB.y, y_max + 0.001f);
        node->BB.z = glm::max(node->BB.z, z_max + 0.001f);
    }
    // if number of triangles < num_node, arrive at the leaf node
    if ((right_index - left_index + 1) <= leaf_num)
    {
        node->num = right_index - left_index + 1;
        node->index = left_index;
        return node;
    }
    float Cost = INF;
    int Axis = 0;
    int Split = (left_index + right_index) / 2;
    for (int axis = 0; axis < 3; axis++)
    {
        // sort according to x/y/z axis
        if (axis == 0)
            std::sort(triangles.begin() + left_index, triangles.begin() + right_index + 1, cmpx);
        if (axis == 1)
            std::sort(triangles.begin() + left_index, triangles.begin() + right_index + 1, cmpy);
        if (axis == 2)
            std::sort(triangles.begin() + left_index, triangles.begin() + right_index + 1, cmpz);

        std::vector<vec3> left_max_point(right_index - left_index + 1, vec3(-INF, -INF, -INF));
        std::vector<vec3> left_min_point(right_index - left_index + 1, vec3(INF, INF, INF));

        for (int i = left_index; i <= right_index; i++)
        {
            Triangle t = triangles[i];
            int bias = (i == left_index) ? 0 : 1;

            left_max_point[i - left_index].x = glm::max(left_max_point[i - left_index - bias].x, glm::max(t.v[0].x, glm::max(t.v[1].x, t.v[2].x)));
            left_max_point[i - left_index].y = glm::max(left_max_point[i - left_index - bias].y, glm::max(t.v[0].y, glm::max(t.v[1].y, t.v[2].y)));
            left_max_point[i - left_index].z = glm::max(left_max_point[i - left_index - bias].z, glm::max(t.v[0].z, glm::max(t.v[1].z, t.v[2].z)));

            left_min_point[i - left_index].x = glm::min(left_min_point[i - left_index - bias].x, glm::min(t.v[0].x, glm::min(t.v[1].x, t.v[2].x)));
            left_min_point[i - left_index].y = glm::min(left_min_point[i - left_index - bias].y, glm::min(t.v[0].y, glm::min(t.v[1].y, t.v[2].y)));
            left_min_point[i - left_index].z = glm::min(left_min_point[i - left_index - bias].z, glm::min(t.v[0].z, glm::min(t.v[1].z, t.v[2].z)));
        }

        std::vector<vec3> right_max_point(right_index - left_index + 1, vec3(-INF, -INF, -INF));
        std::vector<vec3> right_min_point(right_index - left_index + 1, vec3(INF, INF, INF));

        for (int i = right_index; i >= left_index; i--)
        {
            Triangle t = triangles[i];
            int bias = (i == right_index) ? 0 : 1;

            right_max_point[i - left_index].x = glm::max(right_max_point[i - left_index + bias].x, glm::max(t.v[0].x, glm::max(t.v[1].x, t.v[2].x)));
            right_max_point[i - left_index].y = glm::max(right_max_point[i - left_index + bias].y, glm::max(t.v[0].y, glm::max(t.v[1].y, t.v[2].y)));
            right_max_point[i - left_index].z = glm::max(right_max_point[i - left_index + bias].z, glm::max(t.v[0].z, glm::max(t.v[1].z, t.v[2].z)));

            right_min_point[i - left_index].x = glm::min(right_min_point[i - left_index + bias].x, glm::min(t.v[0].x, glm::min(t.v[1].x, t.v[2].x)));
            right_min_point[i - left_index].y = glm::min(right_min_point[i - left_index + bias].y, glm::min(t.v[0].y, glm::min(t.v[1].y, t.v[2].y)));
            right_min_point[i - left_index].z = glm::min(right_min_point[i - left_index + bias].z, glm::min(t.v[0].z, glm::min(t.v[1].z, t.v[2].z)));
        }
        // caculate cost
        float cost = INF;
        int split = left_index;
        for (int i = left_index; i <= right_index - 1; i++)
        {
            float x_length, y_length, z_length;
            vec3 leftAA = left_min_point[i - left_index];
            vec3 leftBB = left_max_point[i - left_index];
            x_length = leftBB.x - leftAA.x;
            y_length = leftBB.y - leftAA.y;
            z_length = leftBB.z - leftAA.z;
            float left_surface_area = 2.0 * ((x_length * y_length) + (x_length * z_length) + (y_length * z_length));
            float leftCost = left_surface_area * (i - left_index + 1);

            vec3 rightAA = right_min_point[i + 1 - left_index];
            vec3 rightBB = right_max_point[i + 1 - left_index];
            x_length = rightBB.x - rightAA.x;
            y_length = rightBB.y - rightAA.y;
            z_length = rightBB.z - rightAA.z;
            float right_surface_area = 2.0 * ((x_length * y_length) + (x_length * z_length) + (y_length * z_length));
            float rightCost = right_surface_area * (right_index - i);

            float totalCost = leftCost + rightCost;
            if (totalCost < cost)
            {
                cost = totalCost;
                split = i;
            }
        }

        if (cost < Cost)
        {
            Cost = cost;
            Axis = axis;
            Split = split;
        }
    }
    // build BVH according to best Axis
    if (Axis == 0)
        std::sort(triangles.begin() + left_index, triangles.begin() + right_index + 1, cmpx);
    if (Axis == 1)
        std::sort(triangles.begin() + left_index, triangles.begin() + right_index + 1, cmpy);
    if (Axis == 2)
        std::sort(triangles.begin() + left_index, triangles.begin() + right_index + 1, cmpz);
    // build left node & right node
    node->left = buildBVH(triangles, left_index, Split, leaf_num);
    node->right = buildBVH(triangles, Split + 1, right_index, leaf_num);

    return node;
}

HitRecord traverseBVH(Ray ray, std::vector<Triangle>& triangles, BVHNode* root)
{
    if (root == NULL)
        return HitRecord();

    if (root->num > 0)
    {
        return interactBVHNode(ray, triangles, root->index, root->index + root->num - 1);
    }

    float d1 = INF, d2 = INF;
    if (root->left)
        d1 = interactAABB(ray, root->left->AA, root->left->BB);
    if (root->right)
        d2 = interactAABB(ray, root->right->AA, root->right->BB);

    HitRecord r1, r2;
    if (d1 > 0)
        r1 = traverseBVH(ray, triangles, root->left);
    if (d2 > 0)
        r2 = traverseBVH(ray, triangles, root->right);
    // if r1.d == r2.d, return light
    if (r1.distance == r2.distance)
        if (r1.triangle.is_emissive)
            return r1;
        else
            return r2;

    return r1.distance < r2.distance ? r1 : r2;
}

HitRecord interactTriangle(Triangle triangle, Ray ray)
{
    HitRecord res;
    vec3 p1 = triangle.v[0], p2 = triangle.v[1], p3 = triangle.v[2];
    vec3 S = ray.startpoint;
    vec3 d = ray.direction;
    vec3 N = triangle.normal;

    if (fabs(dot(N, d)) < 0.00001f)
        return res;
    // caculate distance
    float t = (dot(p1 - S, N)) / dot(d, N);
    if (t < 0.0005f)
        return res;     // if triangle is behing the camera
    vec3 P = S + d * t; // get the hitpoint
    // judege whether the hitpoint is in the triangle
    vec3 c1 = cross(p2 - p1, P - p1);
    vec3 c2 = cross(p3 - p2, P - p2);
    vec3 c3 = cross(p1 - p3, P - p3);
    double dir1 = dot(c1, N), dir2 = dot(c2, N), dir3 = dot(c3, N);
    bool r1 = dir1 > 0 && dir2 > 0 && dir3 > 0;
    bool r2 = dir1 < 0 && dir2 < 0 && dir3 < 0;
    // if hit
    if (r1 || r2)
    {
        res.is_hit = true;
        res.hitpoint = P;
        res.distance = t;
        res.direction = d;
    }

    return res;
}

HitRecord interactBVHNode(Ray ray, std::vector<Triangle> &triangles, int left_index, int right_index)
{
    HitRecord res;
    for (int i = left_index; i <= right_index; i++)
    {
        HitRecord tmp_record = interactTriangle(triangles[i], ray);
        if (tmp_record.is_hit)
        {
            if ((tmp_record.distance == res.distance && triangles[i].is_emissive) || (tmp_record.distance < res.distance))
            { // get the nearest hitpoint int the BVH node
                res = tmp_record;
                res.triangle = triangles[i];
                vec3 barycenter = res.triangle.findBaryCor(res.hitpoint);
                res.pn = normalize((res.triangle.vn[0] * barycenter.x) + (res.triangle.vn[1] * barycenter.y) + (res.triangle.vn[2] * barycenter.z));
            }
        }
    }
    return res;
}

float interactAABB(Ray ray, vec3 AA, vec3 BB)
{
    vec3 inverse_direction = vec3(1.0 / ray.direction.x, 1.0 / ray.direction.y, 1.0 / ray.direction.z);

    vec3 in = (BB - ray.startpoint) * inverse_direction; // ToDo change to /direction
    vec3 out = (AA - ray.startpoint) * inverse_direction;

    vec3 tmax = glm::max(in, out);
    vec3 tmin = glm::min(in, out);

    float t1 = glm::min(tmax.x, glm::min(tmax.y, tmax.z));
    float t0 = glm::max(tmin.x, glm::max(tmin.y, tmin.z));

    return (t1 >= t0) ? ((t0 > 0.0) ? (t0) : (t1)) : (-1);
}

