#include "camera.h"

void Camera::setCamera()
{
    double theta = radians(fovy);
    double h = tan(theta / 2);
    float viewport_height = 2.0 * h;
    float viewport_width = aspect_ratio * viewport_height;

    vec3 w = normalize(eye - lookat);
    vec3 u = normalize(cross(up, w));
    vec3 v = cross(w, u);

    horizontal = viewport_width * u;
    vertical = viewport_height * v;
    lower_left_corner = eye - horizontal / 2.0f - vertical / 2.0f - w;
}

Ray Camera::getRay(float s, float t)
{

    Ray ray;
    ray.startpoint = eye;
    ray.direction = lower_left_corner + s * horizontal + t * vertical - eye;
    ray.direction = normalize(ray.direction);

    return ray;
}

void Camera::Print() {
    printf("Camera:\n");
    printf("fovy: %f ", fovy);
    printf("eye: (%f, %f, %f) ", eye.x, eye.y, eye.z);
    printf("lookat: (%f, %f, %f) ", lookat.x, lookat.y, lookat.z);
    printf("up: (%f, %f, %f) \n", up.x, up.y, up.z);
}