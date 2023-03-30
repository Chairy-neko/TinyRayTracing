#include "pathtracing.h"

vec3 shade(HitRecord &rec, vec3 wi, Scene &scene, BVHNode *root)
{
    vec3 L_dir = vec3(0), L_indir = vec3(0);
    ; // direct light & indirect light

    // for light source
    if (rec.triangle.is_emissive)
    {
        return scene.materials[rec.triangle.mtl_name].radiance; // return light radiance
    }

    // calculate the material infomation
    Material m = scene.materials[rec.triangle.mtl_name];
    vec3 Kd = vec3(0);
    if (m.map_Kd != "")
    { // texture
        vec3 bc = rec.triangle.findBaryCor(rec.hitpoint);
        double col = rec.triangle.vt[0].x * bc.x + rec.triangle.vt[1].x * bc.y + rec.triangle.vt[2].x * bc.z;
        double row = rec.triangle.vt[0].y * bc.x + rec.triangle.vt[1].y * bc.y + rec.triangle.vt[2].y * bc.z;
        double irow = row - floor(row), icol = col - floor(col);
        int r = irow * m.map_height, c = icol * m.map_width;
        cv::Vec3b vec_3 = m.img.at<cv::Vec3b>(r, c);
        Kd.x = (double)vec_3[2] / 255, Kd.y = (double)vec_3[1] / 255, Kd.z = (double)vec_3[0] / 255;
    }
    else
    {
        Kd = m.Kd;
    }

    // direct illumination
    static std::default_random_engine e(time(NULL));
    for (auto light : scene.lights)
    {
        vector<Triangle> light_triangles = scene.materials[light.mtl_name].triangles;
        double total_area = scene.materials[light.mtl_name].area;
        static std::uniform_real_distribution<double> u1(0, total_area);
        double rnd = u1(e);
        for (auto light_triangle : light_triangles)
        {
            if (rnd < light_triangle.area)
            {
                static std::uniform_real_distribution<double> u2(0, 1);
                double rnd1 = u2(e), rnd2 = u2(e), rnd3 = u2(e);
                float p1 = rnd1 / (rnd1 + rnd2 + rnd3), p2 = rnd2 / (rnd1 + rnd2 + rnd3), p3 = rnd3 / (rnd1 + rnd2 + rnd3);
                vec3 light_p = light_triangle.v[0] * p1 + light_triangle.v[1] * p2 + light_triangle.v[2] * p3;
                vec3 light_n = normalize(light_triangle.vn[0] * p1 + light_triangle.vn[1] * p2 + light_triangle.vn[2] * p3);
                vec3 wo = normalize(light_p - rec.hitpoint); // direction from point p to light
                bool visibility = true;                      // judge visibility
                Ray ray_light;
                ray_light.startpoint = rec.hitpoint;
                ray_light.direction = wo;
                HitRecord rec_sample = traverseBVH(ray_light, scene.triangles, root);
                if (rec_sample.triangle.mtl_name != light_triangle.mtl_name)
                {
                    visibility = false;
                }

                if (visibility && dot(wo, rec.pn) > 0)
                {
                    float pdf_light = double(1) / total_area; // pdf_light = 1/A
                    float cos_theta_p = abs(dot(wo, light_n));
                    float cos_theta = abs(dot(wo, rec.pn) / length(rec.pn));
                    vec3 radiance = scene.materials[light_triangle.mtl_name].radiance;
                    vec3 intensity = radiance * cos_theta_p * cos_theta / length2(light_p - rec.hitpoint) / pdf_light;
                    vec3 h = normalize((wi + wo) * 0.5f);
                    double cos_alpha = fmax(dot(rec.pn, h), 0);
                    L_dir += intensity * (Kd / PI + m.Ks * (m.Ns + 2.0f) * (float)pow(cos_alpha, m.Ns) / (2.0f * PI));
                }
                break;
            }
        }
    }

    // return L_dir;
    // indirect illumination
    if (RR(P_RR))
    {
        Ray r = nextRay(rec, -wi, scene);                      // sample next ray
        HitRecord ret = traverseBVH(r, scene.triangles, root); // find the first object ray hits
        if (ret.is_hit && r.ray_type != INVALID)
        {
            vec3 intensity = shade(ret, -r.direction, scene, root) / P_RR;
            switch (r.ray_type)
            {
            case DIFFUSE:
                if (!ret.triangle.is_emissive)
                    L_indir += Kd * intensity;
                break;
            case SPECULAR:
                if (!ret.triangle.is_emissive)
                    L_indir += Kd * intensity;
                break;
            default:
                L_indir += m.Tr * intensity;
            }
        }
    }

    return L_dir + L_indir;
}

bool RR(double prr)
{
    static std::default_random_engine e;
    static std::uniform_real_distribution<double> u1(0, 1);
    return u1(e) < prr;
}

vec3 Sample(vec3 &direction, int ray_type, double Ns)
{
    static std::default_random_engine e(time(NULL));
    static std::uniform_real_distribution<double> u1(0, 1);

    double phi = u1(e) * 2 * PI;
    double theta;
    if (ray_type == DIFFUSE)
    {
        theta = asin(sqrt(u1(e)));
    }
    else if (ray_type == SPECULAR)
    {
        theta = acos(pow(u1(e), (double)1 / (Ns + 1)));
    }
    else
    {
        cerr << "unknown sample ray_type" << endl;
        exit(1);
    }
    vec3 sample(sin(theta) * cos(phi), cos(theta), sin(theta) * sin(phi));
    vec3 front = vec3(0);
    if (fabs(direction.x) > fabs(direction.y))
    {
        front = vec3(direction.z, 0, -direction.x);
        front = normalize(front);
    }
    else
    {
        front = vec3(0, -direction.z, direction.y);
        front = normalize(front);
    }
    vec3 right = cross(direction, front);
    return normalize((right * sample.x) + (direction * sample.y) + (front * sample.z)); // sample_direction;
}

Ray nextRay(HitRecord &rec, vec3 ray_direction, Scene &scene)
{
    static std::default_random_engine e(time(NULL));
    static std::uniform_real_distribution<double> u2(0, 1);

    Material m = scene.materials[rec.triangle.mtl_name];
    vec3 next_ray_direction;

    if (m.Ni > 1)
    {
        double n1, n2;
        double cos_in = dot(ray_direction, rec.pn);
        vec3 normal;
        if (cos_in > 0)
        {
            normal = -rec.pn;
            n1 = m.Ni;
            n2 = 1.0;
        }
        else
        {
            normal = rec.pn;
            n1 = 1.0;
            n2 = m.Ni;
        }

        double rf0 = pow((n1 - n2) / (n1 + n2), 2);
        double fresnel = rf0 + (1.0f - rf0) * pow(1.0f - std::abs(cos_in), 5);
        if (fresnel < u2(e))
        {
            next_ray_direction = refract(ray_direction, normal, (float)(n1 / n2));
            if (next_ray_direction != vec3(0))
            {
                return Ray(rec.hitpoint, next_ray_direction, TRANSMISSION);
            }
            else
            {
                vec3 reflect_ray = reflect(ray_direction, normal); // reflect ray next_ray_direction
                m.Ks = m.Kd;
                return Ray(rec.hitpoint, reflect_ray, SPECULAR);
            }
        }
    }

    double Kd_len = length(m.Kd), Ks_len = length(m.Ks);
    double kd = Kd_len / (Kd_len + Ks_len), ks = Ks_len / (Kd_len + Ks_len);
    int ray_type = INVALID;
    double p = u2(e);
    if (p < kd)
    {
        next_ray_direction = Sample(rec.pn, DIFFUSE, m.Ns);
        ray_type = DIFFUSE;
    }
    else if (m.Ns > 1 && p < kd + ks)
    {
        vec3 reflect_ray = reflect(ray_direction, rec.pn); // reflect ray next_ray_direction
        next_ray_direction = Sample(reflect_ray, SPECULAR, m.Ns);
        ray_type = SPECULAR;
    }

    Ray next_ray(rec.hitpoint, next_ray_direction, ray_type);
    return next_ray;
}