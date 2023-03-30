#include "triangle.h"

double Triangle::calAera()
{
    double a = length(p[1] - p[0]), b = length(p[2] - p[0]), c = length(p[2] - p[1]); // length of a,b,c
    double cos_c = (a * a + b * b - c * c) / (2 * a * b);                             // ”‡œ“∂®¿Ì
    double sin_c = sqrt(1 - pow(cos_c, 2));
    double aera = a * b * sin_c / 2;
    return aera;
}

vec3 Triangle::findBaryCor(vec3 hitp)
{
    Eigen::Matrix<double, 4, 3> A;
    Eigen::Matrix<double, 4, 1> B;
    Eigen::MatrixXd res;

    A << p[0].x, p[1].x, p[2].x,
        p[0].y, p[1].y, p[2].y,
        p[0].z, p[1].z, p[2].z,
        1, 1, 1;
    B << hitp.x, hitp.y, hitp.z, 1;

    res = A.colPivHouseholderQr().solve(B);
    vec3 ret;
    ret.x = res(0, 0), ret.y = res(1, 0), ret.z = res(2, 0);

    return ret;
}