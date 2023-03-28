#ifndef UTILS_H
#define UTILS_H
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <random>

using namespace glm;
using namespace std;

// 0-1 ���������
std::uniform_real_distribution<> dis(0.0, 1.0);
random_device rd;
mt19937 gen(rd());
double randf()
{
    return dis(gen);
}

// ��λ���ڵ��������
vec3 random_in_unit_sphere()
{

    vec3 d;
    do
    {
        d = 2.0f * vec3(randf(), randf(), randf()) - vec3(1, 1, 1);
    } while (dot(d, d) > 1.0);
    return normalize(d);
}

// ��������������
vec3 random_in_hemisphere(vec3 n)
{
    // ������
    return normalize(random_in_unit_sphere() + n);
}

#endif