#pragma once
#include <random>
#include <glm/glm.hpp>
#include <iostream>

using namespace std;
using namespace glm;

const float PI = 3.1415926f;
const float P_RR = 0.8; // russian roulette probability

int num_kd = 0, num_ks = 0; // for debugging purposes
int dark = 0;

bool russian_Roulette(double pr)
{
    static std::default_random_engine e;
    static std::uniform_real_distribution<double> u1(0, 1);
    double rnd = u1(e);
    // std::cout << rnd << std::endl;
    if (rnd < pr)
        return true;
    else
        return false;
}

// importance sampling based on BRDF
vec3 BRDFImportanceSampling(vec3 &direction, int type, double Ns)
{
    // Ref : https://inst.eecs.berkeley.edu/~cs283/sp13/lectures/283-lecture11.pdf
    static std::default_random_engine e(time(NULL));
    static std::uniform_real_distribution<double> u1(0, 1);

    double phi = u1(e) * 2 * PI;
    double theta;
    if (type == DIFFUSE)
    { // diffuse
        // std::cout << "sample diffuse :" << f.material << std::endl;
        //  for cosine-weighted Lambertian
        theta = asin(sqrt(u1(e)));
    }
    else if (type == SPECULAR)
    {
        // for sampling specular term
        // std::cout << "sample speculat :" << f.material << std::endl;
        theta = acos(pow(u1(e), (double)1 / (Ns + 1)));
    }
    else
    {
        cerr << "unknown sample type" << endl;
    }
    vec3 sample(sin(theta) * cos(phi), cos(theta), sin(theta) * sin(phi));
    vec3 front;
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
    vec3 ret = (right * sample.x) + (direction * sample.y) + (front * sample.z);
    ret = normalize(ret);
    return ret;
}

Ray nextRay(HitResult &res, vec3 dir, Scene &scene)
{
    static std::default_random_engine e(time(NULL));
    static std::uniform_real_distribution<double> u2(0, 1);

    Material m = scene.materials[res.triangle.mtl_name];

    vec3 direction;

    // for fraction

    if (m.Ni > 1)
    {
        // std::cout << "fraction :" << m->name << std::endl;
        double n1, n2;
        double cos_in = dot(-dir, res.pn);
        vec3 normal;
        if (cos_in > 0)
        {
            // out of glass
            // std::cout << "out" << std::endl;
            normal = -res.pn;
            n1 = m.Ni;
            n2 = 1.0;
        }
        else
        {
            // std::cout << "in" << std::endl;
            //  in to the glass
            normal = res.pn;
            n1 = 1.0;
            n2 = m.Ni;
        }

        // Ref : https://en.wikipedia.org/wiki/Schlick%27s_approximation
        double rf0 = pow((n1 - n2) / (n1 + n2), 2);
        double fresnel = rf0 + (1.0f - rf0) * pow(1.0f - std::abs(cos_in), 5);
        if (fresnel < u2(e))
        {
            direction = refract(-dir, normal, (float)(n1 / n2));
            if (direction != vec3(0))
            {
                return Ray(res.hitpoint, direction, TRANSMISSION);
            }
            else
            {
                // printf("all specular\n");
                // only specular for refraction material
                // vec3 incoming = -dir;
                vec3 reflect_ray = reflect(-dir, normal); // reflect ray direction
                // reflect = incoming - normal * dot(incoming, normal) * 2.0f;
                return Ray(res.hitpoint, reflect_ray, SPECULAR);
                // only specular for refraction material
            }
        }
    }

    // double Kd_len = length(Kd), Ks_len = length(m.Ks);
    double kd = (m.Kd.x + m.Kd.y + m.Kd.z) / 3.0f, ks = (m.Ks.x + m.Ks.y + m.Ks.z) / 3.0f;
    int ray_type = INVALID;
    double p = u2(e);
    if (p < kd)
    {
        num_kd++;
        // sample diffuse
        direction = BRDFImportanceSampling(res.pn, DIFFUSE, m.Ns);
        ray_type = DIFFUSE;
    }
    else if (m.Ns > 1 && p < kd + ks)
    {
        num_ks++;
        // sample specular
        //  vec3 incoming = -dir;
        vec3 reflect_ray = reflect(-dir, res.pn); // reflect ray direction
        // reflect_ray = incoming - p.pn * (incoming * p.pn) * 2;
        direction = BRDFImportanceSampling(reflect_ray, SPECULAR, m.Ns);
        ray_type = SPECULAR;
    }

    Ray r(res.hitpoint + (direction * 0.001f), direction, ray_type);
    return r;
}

vec3 shade(HitResult &res,
           vec3 dir,
           Scene &scene,
           BVHNode *root)
{
    vec3 L_dir = vec3(0); // direct light

    if (res.triangle.isEmissive)
    {                                                           // for light
        return scene.materials[res.triangle.mtl_name].radiance; // return light radiance
    }

    // calculate the material infomation
    Material m = scene.materials[res.triangle.mtl_name];
    vec3 Kd = vec3(0), Ks = m.Ks;
    if (m.map_Kd != "")
    { // texture
        vec3 bc = res.triangle.findBaryCor(res.hitpoint);
        // garcov.print();
        // double row = res.triangle.vt[0].x * bc.x + res.triangle.vt[1].x * bc.y + res.triangle.vt[2].x * bc.z; //
        double col = res.triangle.vt[0].x * bc.x + res.triangle.vt[1].x * bc.y + res.triangle.vt[2].x * bc.z; //
        // std::cout << p.f.vt1.vtx << " " << p.f.vt2.vtx << " " << p.f.vt3.vtx << std::endl;
        // double col = res.triangle.vt[0].y * bc.x + res.triangle.vt[1].y * bc.y + res.triangle.vt[2].y * bc.z;
        double row = res.triangle.vt[0].y * bc.x + res.triangle.vt[1].y * bc.y + res.triangle.vt[2].y * bc.z;
        double irow = row - floor(row), icol = col - floor(col);
        int r = irow * m.map_height, c = icol * m.map_width;
        // std::cout << m->name << " " << r << " " << c << std::endl;
        cv::Vec3b vec_3 = m.img.at<cv::Vec3b>(r, c);
        Kd.x = (double)vec_3[2] / 255, Kd.y = (double)vec_3[1] / 255, Kd.z = (double)vec_3[0] / 255;
    }
    else
    {
        Kd = m.Kd;
    }

    // direct illumination
    // Uniformly sample light at x'
    static std::default_random_engine e(time(NULL));
    Triangle sample_face;
    for (int i = 0; i < scene.lights.size(); i++)
    {
        // for every direct light
        vec3 xl = vec3(0);         // xl is the sample point on light
        vec3 vn = vec3(0);         // normal at point xl
        Light l = scene.lights[i]; // l record the light attribute

        double total_aera = scene.materials[l.mtl_name].area;
        vector<Triangle> light_triangles = scene.materials[l.mtl_name].triangles;
        // int n_triangle_face = light_triangles.size();
        // double *triangle_aera = new double[n_triangle_face];
        // for (int j = 0; j < n_triangle_face; j++)
        // {
        //     // for every triangle mesh
        //     total_aera += light_triangles[j].area; // this is for pdf_light
        //     triangle_aera[j] = total_aera;
        // }
        static std::uniform_real_distribution<double> u1(0, total_aera);
        double rnd = u1(e);
        // std::cout << rnd << std::endl;
        // uniformly generate a random number based on aera to choose which triagnle mesh to sample
        for (auto light_triangle : light_triangles)
        {
            if (rnd < light_triangle.area)
            {
                // sample a vertex on triangle mesh j
                sample_face = light_triangle;
                static std::uniform_real_distribution<double> u2(0, 1);
                double rnd1 = u2(e), rnd2 = u2(e), rnd3 = u2(e);
                float p1 = rnd1 / (rnd1 + rnd2 + rnd3), p2 = rnd2 / (rnd1 + rnd2 + rnd3), p3 = rnd3 / (rnd1 + rnd2 + rnd3);
                xl = sample_face.p[0] * p1 + sample_face.p[1] * p2 + sample_face.p[2] * p3;
                vn = normalize(sample_face.vn[0] * p1 + sample_face.vn[1] * p2 + sample_face.vn[2] * p3);
                break;
            }
        }

        vec3 direction = normalize(xl - res.hitpoint); // direction from point p to light
        float visibility = 1;                          // judge visibility
        Ray rl;
        rl.startPoint = res.hitpoint + (direction * 0.001f); // add micro turbulance
        rl.direction = direction;
        HitResult res2 = hitBVH(rl, scene.triangles, root);
        if (res2.triangle.mtl_name != sample_face.mtl_name)
        {
            visibility = 0;
        }

        if (dot(direction, res.pn) > 0)
        {
            float pdf_light = double(1) / total_aera; // pdf_light = 1/A where A is the aera of light source
            float cos_theta = abs(dot(direction, vn));
            // float cos_theta = abs(dot(direction, vn) / length(direction) / length(vn));
            float cos_theta_hat = abs(dot(direction, res.pn) / length(res.pn));
            // float cos_theta_hat = abs(dot(direction, res.pn) / length(direction) / length(res.pn));
            vec3 intensity = l.radiance * cos_theta * cos_theta_hat / length2(xl - res.hitpoint) / pdf_light * visibility;
            // float kd_dots = dot(direction, res.pn); // cos between light to intersection and face normal
            vec3 reflect_ray = reflect(direction, res.pn);
            double cos_alpha = fmax(dot(dir, reflect_ray), 0.0f);
            L_dir += intensity * (Kd / PI + Ks * (m.Ns + 2.0f) * (float)pow(cos_alpha, m.Ns) / (2.0f * PI));
            // only add diffuse
            // if (kd_dots > 0)
            // {
            //     L_dir += Kd * intensity * kd_dots / PI;
            //     // L_dir.x += Kd.x * intensity.x * kd_dots / PI;
            //     // L_dir.y += Kd.y * intensity.y * kd_dots / PI;
            //     // L_dir.z += Kd.z * intensity.z * kd_dots / PI;
            // }
        }
    }

    // return L_dir;
    //   indirect illumination
    vec3 L_indir = vec3(0);

    // BRDF sample the hemisphere, get direction
    if (russian_Roulette(P_RR))
    {
        Ray r = nextRay(res, dir, scene); // sample next ray based on BRDF
        HitResult ret = hitBVH(r, scene.triangles, root);// find the first object ray hits
        if (ret.isHit && r.ray_type != INVALID)
        {
            vec3 intensity = shade(ret, -r.direction, scene, root) / P_RR;
            if (r.ray_type == DIFFUSE)
            {
                if (!res.triangle.isEmissive)
                { // hit none emitting object
                    L_indir += Kd * intensity;
                    // L_indir.x += kd.x * intensity.x;
                    // L_indir.y += kd.y * intensity.y;
                    // L_indir.z += kd.z * intensity.z;
                }
            }
            else if (r.ray_type == SPECULAR)
            {
                if (!res.triangle.isEmissive)
                {
                    L_indir += m.Ks * intensity;
                    // L_indir.x += m->ks.x *  intensity.x;
                    // L_indir.y += m->ks.y *  intensity.y;
                    // L_indir.z += m->ks.z *  intensity.z;
                }
            }
            else
            {
                //printf("TRANSMISSTION\n");
                L_indir += m.Tr * intensity;
            }
        }
    }
    if (L_indir.x < 0 || L_indir.y < 0 || L_indir.z < 0)
        printf("deubg dark: %d\n", dark++);

    return L_dir + L_indir;
}
