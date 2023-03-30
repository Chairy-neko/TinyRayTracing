#include "triangle.h"

double Triangle::calAera()
{
    double a = length(v[1] - v[0]), b = length(v[2] - v[0]), c = length(v[2] - v[1]);
    double cos_c = (a * a + b * b - c * c) / (2 * a * b);
    double sin_c = sqrt(1 - pow(cos_c, 2));
    double aera = a * b * sin_c / 2;
    return aera;
}

vec3 Triangle::findBaryCor(vec3 hitp)
{
    Eigen::Matrix<double, 4, 3> A;
    Eigen::Matrix<double, 4, 1> B;
    Eigen::MatrixXd res;

    A << v[0].x, v[1].x, v[2].x,
        v[0].y, v[1].y, v[2].y,
        v[0].z, v[1].z, v[2].z,
        1, 1, 1;
    B << hitp.x, hitp.y, hitp.z, 1;

    res = A.colPivHouseholderQr().solve(B);
    vec3 ret;
    ret.x = res(0, 0), ret.y = res(1, 0), ret.z = res(2, 0);

    return ret;
}