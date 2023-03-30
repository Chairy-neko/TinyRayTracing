#include "pathtracing.h"

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
                return Ray(res.hitpoint, reflect_ray, TREFLECTION);
                // only specular for refraction material
            }
        }
    }

    double Kd_len = length(m.Kd), Ks_len = length(m.Ks);
    double kd = Kd_len / (Kd_len + Ks_len), ks = Ks_len / (Kd_len + Ks_len);
    // double kd = (m.Kd.x + m.Kd.y + m.Kd.z) / 3.0f, ks = (m.Ks.x + m.Ks.y + m.Ks.z) / 3.0f;
    int ray_type = INVALID;
    double p = u2(e);
    if (p < kd)
    {
        // num_kd++;
        // sample diffuse
        direction = BRDFImportanceSampling(res.pn, DIFFUSE, m.Ns);
        ray_type = DIFFUSE;
    }
    else if (m.Ns > 1 && p < kd + ks)
    {
        // num_ks++;
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

vec3 directIllumination(HitResult &res, vec3 dir, Scene &scene, BVHNode *root, vec3 Kd)
{
    Material m = scene.materials[res.triangle.mtl_name];
    vec3 L_dir = vec3(0);
    static std::default_random_engine e(time(NULL));
    // for (auto light : scene.lights)
    // {
    //     vec3 light_p = vec3(0);
    //     vec3 light_n = vec3(0);

    //     double total_aera = scene.materials[light.mtl_name].area;

    //     vector<Triangle> light_triangles = scene.materials[light.mtl_name].triangles;
    //     static std::uniform_real_distribution<double> u1(0, total_aera);
    //     double rnd = u1(e);
    //     for (auto light_triangle : light_triangles)
    //     {
    //         if (rnd < light_triangle.area)
    //         {
    //             static std::uniform_real_distribution<double> u2(0, 1);
    //             double rnd1 = u2(e), rnd2 = u2(e), rnd3 = u2(e);
    //             float p1 = rnd1 / (rnd1 + rnd2 + rnd3), p2 = rnd2 / (rnd1 + rnd2 + rnd3), p3 = rnd3 / (rnd1 + rnd2 + rnd3);
    //             light_p = light_triangle.p[0] * p1 + light_triangle.p[1] * p2 + light_triangle.p[2] * p3;
    //             light_n = normalize(light_triangle.vn[0] * p1 + light_triangle.vn[1] * p2 + light_triangle.vn[2] * p3);
    //             vec3 direction = normalize(light_p - res.hitpoint); // direction from point p to light
    //             float visibility = 1;                               // judge visibility
    //             Ray rl;
    //             rl.startPoint = res.hitpoint + (direction * 0.001f); // add micro turbulance
    //             rl.direction = direction;
    //             HitResult res2 = hitBVH(rl, scene.triangles, root);
    //             if (res2.triangle.mtl_name != light_triangle.mtl_name)
    //             {
    //                 visibility = 0;
    //             }

    //             if (dot(direction, res.pn) > 0)
    //             {
    //                 float pdf_light = double(1) / scene.materials[light_triangle.mtl_name].area; // pdf_light = 1/A where A is the aera of light source
    //                 float cos_theta = abs(dot(direction, light_n));
    //                 // float cos_theta = abs(dot(direction, light_n) / length(direction) / length(light_n));
    //                 float cos_theta_hat = abs(dot(direction, res.pn) / length(res.pn));
    //                 // float cos_theta_hat = abs(dot(direction, res.pn) / length(direction) / length(res.pn));
    //                 vec3 intensity = light.radiance * cos_theta * cos_theta_hat / length2(light_p - res.hitpoint) / pdf_light * visibility;
    //                 // float kd_dots = dot(direction, res.pn); // cos between light to intersection and face normal
    //                 L_dir += intensity * Kd / PI;
    //                 vec3 h = normalize((dir + direction) * 0.5f);
    //                 double cos_alpha = dot(res.pn, h);
    //                 // vec3 reflect_ray = reflect(-direction, res.pn);
    //                 // vec3 reflect_ray = reflect(dir, res.pn);
    //                 // double cos_alpha = fmax(dot(dir, reflect_ray), 0.0f);
    //                 // L_dir += intensity * (Kd / PI);
    //                 if (dot(res.pn, direction) > 0)
    //                     L_dir += intensity * (m.Ks * (m.Ns + 2.0f) * (float)pow(cos_alpha, m.Ns) / (2.0f * PI));
    //             }
    //             break;
    //         }
    //     }
    // }

    // area
    // double total_area = scene.total_area;
    // static std::uniform_real_distribution<double> u1(0, total_area);
    // // for (auto light : scene.lights)
    // // {
    // //     vector<Triangle> light_triangle_tmp = scene.materials[light.mtl_name].triangles;
    // //     light_triangles.insert(light_triangles.end(), light_triangle_tmp.begin(), light_triangle_tmp.end());
    // // }
    // bool flag = false;
    // double rnd = u1(e);
    // for (auto light : scene.lights)
    // {
    //     vector<Triangle> light_triangles = scene.materials[light.mtl_name].triangles;
    //     for (auto light_triangle : light_triangles)
    //     {
    //         if (rnd < light_triangle.area)
    //         {
    //             static std::uniform_real_distribution<double> u2(0, 1);
    //             double rnd1 = u2(e), rnd2 = u2(e), rnd3 = u2(e);
    //             float p1 = rnd1 / (rnd1 + rnd2 + rnd3), p2 = rnd2 / (rnd1 + rnd2 + rnd3), p3 = rnd3 / (rnd1 + rnd2 + rnd3);
    //             vec3 light_p = light_triangle.p[0] * p1 + light_triangle.p[1] * p2 + light_triangle.p[2] * p3;
    //             vec3 light_n = normalize(light_triangle.vn[0] * p1 + light_triangle.vn[1] * p2 + light_triangle.vn[2] * p3);
    //             vec3 direction = normalize(light_p - res.hitpoint); // direction from point p to light
    //             float visibility = 1;                               // judge visibility
    //             Ray rl;
    //             rl.startPoint = res.hitpoint + (direction * 0.001f); // add micro turbulance
    //             rl.direction = direction;
    //             HitResult res2 = hitBVH(rl, scene.triangles, root);
    //             if (res2.triangle.mtl_name != light_triangle.mtl_name)
    //             {
    //                 visibility = 0;
    //             }

    //             if (dot(direction, res.pn) > 0)
    //             {
    //                 float pdf_light = double(1) / scene.materials[light_triangle.mtl_name].area; // pdf_light = 1/A
    //                 float cos_theta = abs(dot(direction, light_n));
    //                 float cos_theta_hat = abs(dot(direction, res.pn) / length(res.pn));
    //                 vec3 radiance = scene.materials[light_triangle.mtl_name].radiance;
    //                 vec3 intensity = radiance * cos_theta * cos_theta_hat / length2(light_p - res.hitpoint) / pdf_light * visibility;
    //                 L_dir += intensity * Kd / PI;
    //                 vec3 h = normalize((dir + direction) * 0.5f);
    //                 double cos_alpha = dot(res.pn, h);
    //                 L_dir += intensity * (m.Ks * (m.Ns + 2.0f) * (float)pow(cos_alpha, m.Ns) / (2.0f * PI));
    //             }
    //             flag = true;
    //             break;
    //         }
    //     }
    //     if (flag)
    //         break;
    // }

    // double total_area_radiance = 0.0;
    // double total_area = 0.0;
    // vector<double> area_radiance;
    // for (auto light : scene.lights)
    // {
    //     Material light_m = scene.materials[light.mtl_name];
    //     total_area += light_m.area;
    //     total_area_radiance += length(light_m.radiance);
    //     area_radiance.push_back(total_area_radiance);
    // }
    // static std::uniform_real_distribution<double> u1(0, scene.total_area);
    // double rnd_mtl = u1(e);
    // for (int i = 0; i < scene.lights.size(); i++)
    // {
    //     if (rnd_mtl < scene.materials[scene.lights[i].mtl_name].area)
    //     {
    //         Light light = scene.lights[i];
    //         static std::uniform_real_distribution<double> u2(0, scene.materials[light.mtl_name].area);
    //         double rnd_tri = u2(e);
    //         vector<Triangle> light_triangles = scene.materials[light.mtl_name].triangles;
    //         for (auto light_triangle : light_triangles)
    //         {
    //             vec3 radiance = scene.materials[light_triangle.mtl_name].radiance;
    //             if (rnd_tri < light_triangle.area)
    //             {
    //                 static std::uniform_real_distribution<double> u2(0, 1);
    //                 double rnd1 = u2(e), rnd2 = u2(e), rnd3 = u2(e);
    //                 float p1 = rnd1 / (rnd1 + rnd2 + rnd3), p2 = rnd2 / (rnd1 + rnd2 + rnd3), p3 = rnd3 / (rnd1 + rnd2 + rnd3);
    //                 vec3 light_p = light_triangle.p[0] * p1 + light_triangle.p[1] * p2 + light_triangle.p[2] * p3;
    //                 vec3 light_n = normalize(light_triangle.vn[0] * p1 + light_triangle.vn[1] * p2 + light_triangle.vn[2] * p3);
    //                 vec3 direction = normalize(light_p - res.hitpoint); // direction from point p to light
    //                 float visibility = 1;                               // judge visibility
    //                 Ray rl;
    //                 rl.startPoint = res.hitpoint + (direction * 0.001f); // add micro turbulance
    //                 rl.direction = direction;
    //                 HitResult res2 = hitBVH(rl, scene.triangles, root);
    //                 if (res2.triangle.mtl_name != light_triangle.mtl_name)
    //                 {
    //                     visibility = 0;
    //                 }

    //                 float pdf_light = double(1) / total_area; // pdf_light = 1/A
    //                 float cos_theta = abs(dot(direction, light_n));
    //                 float cos_theta_hat = abs(dot(direction, res.pn));
    //                 vec3 intensity = radiance * cos_theta * cos_theta_hat / length2(light_p - res.hitpoint) / pdf_light * visibility;
    //                 L_dir += intensity * Kd / PI;
    //                 if (dot(direction, res.pn) > 0)
    //                 {
    //                     vec3 h = normalize((dir + direction) * 0.5f);
    //                     double cos_alpha = dot(res.pn, h);
    //                     L_dir += intensity * (m.Ks * (m.Ns + 2.0f) * (float)pow(cos_alpha, m.Ns) / (2.0f * PI));
    //                 }
    //                 break;
    //             }
    //         }
    //         break;
    //     }
    // }

    int light_triangle_num = scene.light_triangles.size();
    static std::uniform_real_distribution<double> u1(0, light_triangle_num);
    int rnd = floor(u1(e));
    Triangle light_triangle = scene.light_triangles[rnd];
    static std::uniform_real_distribution<double> u2(0, 1);
    double rnd1 = u2(e), rnd2 = u2(e), rnd3 = u2(e);
    float p1 = rnd1 / (rnd1 + rnd2 + rnd3), p2 = rnd2 / (rnd1 + rnd2 + rnd3), p3 = rnd3 / (rnd1 + rnd2 + rnd3);
    vec3 light_p = light_triangle.p[0] * p1 + light_triangle.p[1] * p2 + light_triangle.p[2] * p3;
    vec3 light_n = normalize(light_triangle.vn[0] * p1 + light_triangle.vn[1] * p2 + light_triangle.vn[2] * p3);
    vec3 direction = normalize(light_p - res.hitpoint); // direction from point p to light
    float visibility = 1;                               // judge visibility
    Ray rl;
    rl.startPoint = res.hitpoint + (direction * 0.001f); // add micro turbulance
    rl.direction = direction;
    HitResult res2 = hitBVH(rl, scene.triangles, root);
    if (res2.triangle.mtl_name != light_triangle.mtl_name)
    {
        visibility = 0;
    }

    if (dot(direction, res.pn) > 0)
    {
        float pdf_light = double(1) / scene.materials[light_triangle.mtl_name].area; // pdf_light = 1/A
        float cos_theta = abs(dot(direction, light_n));
        float cos_theta_hat = abs(dot(direction, res.pn) / length(res.pn));
        vec3 radiance = scene.materials[light_triangle.mtl_name].radiance;
        vec3 intensity = radiance * cos_theta * cos_theta_hat / length2(light_p - res.hitpoint) / pdf_light * visibility;
        L_dir += intensity * Kd / PI;
        vec3 h = normalize((dir + direction) * 0.5f);
        double cos_alpha = dot(res.pn, h);
        L_dir += intensity * (m.Ks * (m.Ns + 2.0f) * (float)pow(cos_alpha, m.Ns) / (2.0f * PI));
    }

    return L_dir;
}

vec3 shade(HitResult &res, vec3 dir, Scene &scene, BVHNode *root)
{
    vec3 L_dir = vec3(0); // direct light

    if (res.triangle.isEmissive)
    {                                                           // for light
        return scene.materials[res.triangle.mtl_name].radiance; // return light radiance
    }

    // calculate the material infomation
    Material m = scene.materials[res.triangle.mtl_name];
    vec3 Kd = vec3(0);
    if (m.map_Kd != "")
    { // texture
        vec3 bc = res.triangle.findBaryCor(res.hitpoint);
        double col = res.triangle.vt[0].x * bc.x + res.triangle.vt[1].x * bc.y + res.triangle.vt[2].x * bc.z;
        double row = res.triangle.vt[0].y * bc.x + res.triangle.vt[1].y * bc.y + res.triangle.vt[2].y * bc.z;
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
    L_dir = directIllumination(res, dir, scene, root, Kd);

    // return L_dir;
    //    indirect illumination
    vec3 L_indir = vec3(0);

    // BRDF sample the hemisphere, get direction
    if (russian_Roulette(P_RR))
    {
        Ray r = nextRay(res, dir, scene);                 // sample next ray based on BRDF
        HitResult ret = hitBVH(r, scene.triangles, root); // find the first object ray hits
        if (ret.isHit && r.ray_type != INVALID)
        {
            vec3 intensity = shade(ret, -r.direction, scene, root) / P_RR;
            if (r.ray_type == DIFFUSE)
            {
                if (!ret.triangle.isEmissive)
                { // hit none emitting object
                    num_kd++;
                    L_indir += Kd * intensity;
                    // L_indir.x += kd.x * intensity.x;
                    // L_indir.y += kd.y * intensity.y;
                    // L_indir.z += kd.z * intensity.z;
                }
            }
            else if (r.ray_type == SPECULAR)
            {
                if (!ret.triangle.isEmissive)
                {
                    num_ks++;
                    L_indir += m.Ks * intensity;
                    // L_indir.x += m->ks.x *  intensity.x;
                    // L_indir.y += m->ks.y *  intensity.y;
                    // L_indir.z += m->ks.z *  intensity.z;
                }
            }
            else if (r.ray_type == TREFLECTION)
            {
                L_indir += m.Kd * intensity;
            }
            else
            {
                // printf("TRANSMISSTION\n");
                L_indir += m.Tr * intensity;
            }
        }
    }
    if (L_indir.x < 0 || L_indir.y < 0 || L_indir.z < 0)
        printf("deubg dark: %d\n", dark++);

    return L_dir + L_indir;
}
