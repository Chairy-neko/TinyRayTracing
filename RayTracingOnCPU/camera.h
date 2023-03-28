#pragma once
#include "utils.h"
#include "ray.h"

//// 相机参数
//const double SCREEN_Z = 1.1;        // 视平面 z 坐标
//const vec3 EYE = vec3(0, 0, 4.0);   // 相机位置

class Camera
{
public:
	Camera() {}
	~Camera() {}
	double fovy = 90;
	vec3 eye = vec3(278.0, 273.0, -800.0);
	vec3 lookat = vec3(278.0, 273.0, -799.0);
	vec3 up = vec3(0.0, 1.0, 0.0);
	double aspect_ratio = 1.0;

	vec3 lower_left_corner = vec3(0.0, 0.0, 0.0);
	vec3 horizontal = vec3(0.0, 0.0, 0.0);
	vec3 vertical = vec3(0.0, 0.0, 0.0);

	void SetCamera(double a) {
		aspect_ratio = a;
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

	Ray get_ray(float s, float t) const{

		Ray ray;
		ray.startPoint = eye;
		ray.direction = lower_left_corner + s * horizontal + t * vertical - eye;
		ray.direction = normalize(ray.direction);

		return ray;
	}
};