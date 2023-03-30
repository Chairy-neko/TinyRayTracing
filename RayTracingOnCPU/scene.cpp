#include "scene.h"

void Scene::readxml(string xml_path) // get camera & lights
{
    XMLDocument doc;
    XMLError error = doc.LoadFile(xml_path.c_str());
    if (error != XML_SUCCESS)
    {
        std::cout << "Read xml failed." << doc.ErrorStr() << std::endl;
        exit(1);
    }
    XMLElement *xml_camera = doc.RootElement();
    img_width = stoi(xml_camera->Attribute("width"));
    img_height = stoi(xml_camera->Attribute("height"));
    camera.aspect_ratio = (double)img_width / (double)img_height;
    camera.fovy = stof(xml_camera->Attribute("fovy"));
    XMLElement *xml_eye = xml_camera->FirstChildElement("eye");
    camera.eye = vec3(stof(xml_eye->Attribute("x")), stof(xml_eye->Attribute("y")), stof(xml_eye->Attribute("z")));
    XMLElement *xml_lookat = xml_camera->FirstChildElement("lookat");
    camera.lookat = vec3(stof(xml_lookat->Attribute("x")), stof(xml_lookat->Attribute("y")), stof(xml_lookat->Attribute("z")));
    XMLElement *xml_up = xml_camera->FirstChildElement("up");
    camera.up = vec3(stof(xml_up->Attribute("x")), stof(xml_up->Attribute("y")), stof(xml_up->Attribute("z")));
    XMLElement *xml_light = xml_camera->NextSiblingElement("light");
    camera.SetCamera();
    while (xml_light)
    {
        string mtl_name = xml_light->Attribute("mtlname");
        string radiance_str = xml_light->Attribute("radiance");
        vec3 radiance = vec3(0.0, 0.0, 0.0);
        int tmp = 0, flag = 0;
        for (int i = 0; i < radiance_str.length(); i++)
        {
            if (radiance_str[i] == ',' && flag == 0)
            {
                radiance.x = stof(radiance_str.substr(tmp, i - tmp));
                tmp = i + 1;
                flag = 1;
            }
            else if (radiance_str[i] == ',' && flag == 1)
            {
                radiance.y = stof(radiance_str.substr(tmp, i - tmp));
                tmp = i + 1;
                flag = 2;
            }
            else if (i == radiance_str.length() - 1)
            {
                radiance.z = stof(radiance_str.substr(tmp, i - tmp + 1));
            }
        }
        lights.push_back(Light(mtl_name, radiance));
        materials[mtl_name].isEmissive = true;
        materials[mtl_name].radiance = radiance;
        xml_light = xml_light->NextSiblingElement("light");
    }
}

void Scene::readmtl(string mtl_path, string basedir)
{
    std::ifstream fin(mtl_path);
    std::string line;
    if (!fin.is_open())
    {
        std::cout << "Read " << mtl_path << " failed." << std::endl;
        exit(-1);
    }
    string mtlname = "";
    while (std::getline(fin, line))
    {
        std::istringstream sin(line);
        std::string type;
        string mapkd = "";
        float x, y, z;
        float n;

        sin >> type;
        if (type == "newmtl")
        {
            sin >> mtlname;
        }
        else if (type == "Kd")
        {
            sin >> x >> y >> z;
            materials[mtlname].Kd = vec3(x, y, z);
        }
        else if (type == "Ks")
        {
            sin >> x >> y >> z;
            materials[mtlname].Ks = vec3(x, y, z);
        }
        else if (type == "Tr")
        {
            sin >> x >> y >> z;
            materials[mtlname].Tr = vec3(x, y, z);
        }
        else if (type == "Ns")
        {
            sin >> n;
            materials[mtlname].Ns = n;
        }
        else if (type == "Ni")
        {
            sin >> n;
            materials[mtlname].Ni = n;
        }
        else if (type == "map_Kd")
        {
            sin >> mapkd;
            materials[mtlname].map_Kd = basedir + "/" + mapkd;
            materials[mtlname].readinMap();
        }
    }
    printf("num of materials: %d\n", (int)materials.size());
}

void Scene::readobj(string obj_path)
{
    std::ifstream fin(obj_path);
    std::string line;
    if (!fin.is_open())
    {
        std::cout << "Read " << obj_path << " failed." << std::endl;
        exit(-1);
    }

    std::vector<vec3> vertices;
    std::vector<vec3> vn;
    std::vector<vec2> vt;
    bool isvnvt = true;
    int triangleid = 0;
    std::string mtl_name = "";
    while (std::getline(fin, line))
    {
        std::istringstream sin(line);
        std::string type;
        float x, y, z;
        string f[3];

        sin >> type;
        if (type == "v")
        {
            sin >> x >> y >> z;
            vertices.push_back(vec3(x, y, z));
        }
        else if (type == "vn")
        {
            sin >> x >> y >> z;
            vn.push_back(vec3(x, y, z));
        }
        else if (type == "vt")
        {
            if (vn.size() == 0)
                isvnvt = false;
            sin >> x >> y;
            vt.push_back(vec2(x, y));
        }
        else if (type == "usemtl")
        {
            sin >> mtl_name;
        }
        else if (type == "f")
        {
            sin >> f[0] >> f[1] >> f[2];
            int tmp = 0, tmpindex = 0;
            Triangle triangle;
            for (int fi = 0; fi < 3; fi++)
            {
                for (int index = 0; index < f[fi].length(); index++)
                {
                    if (f[fi][index] == '/' && tmp == 0)
                    {
                        tmpindex = stoi(f[fi].substr(tmp, index - tmp)) - 1;
                        triangle.p[fi] = vertices[tmpindex];
                        tmp = index + 1;
                    }
                    else if (f[fi][index] == '/' && tmp != 0)
                    {
                        tmpindex = stoi(f[fi].substr(tmp, index - tmp)) - 1;
                        if (isvnvt)
                            triangle.vn[fi] = vn[tmpindex];
                        else
                            triangle.vt[fi] = vt[tmpindex];
                        tmp = index + 1;
                    }
                    else if (index == f[fi].length() - 1)
                    {
                        tmpindex = stoi(f[fi].substr(tmp, index - tmp + 1)) - 1;
                        if (isvnvt)
                            triangle.vt[fi] = vt[tmpindex];
                        else
                            triangle.vn[fi] = vn[tmpindex];
                        tmp = 0;
                    }
                }
                tmp = 0;
            }
            triangle.normal = normalize(cross(triangle.p[1] - triangle.p[0], triangle.p[2] - triangle.p[0]));
            triangle.center = (triangle.p[0] + triangle.p[1] + triangle.p[2]) / vec3(3.0, 3.0, 3.0);
            triangle.mtl_name = mtl_name;
            triangle.id = triangleid++;
            if(materials[mtl_name].isEmissive){
                triangle.isEmissive = true;
                // total_area += triangle.calAera();
                materials[mtl_name].area += triangle.calAera();
                // triangle.area = total_area;
                light_triangles.push_back(triangle);
                // materials[mtl_name].triangles.push_back(triangle);

                // materials[mtl_name].area += triangle.calAera();
                // triangle.area = materials[mtl_name].area;
                // materials[mtl_name].triangles.push_back(triangle);
            }
            triangles.push_back(triangle);
        }
    }
    printf("num of vertices: %d\n", (int)vertices.size());
    printf("num of vn: %d\n", (int)vn.size());
    printf("num of vt: %d\n", (int)vt.size());
    printf("num of triangles: %d\n", (int)triangles.size());
}

//void Camera::SetCamera()
//{
//    double theta = radians(fovy);
//    double h = tan(theta / 2);
//    float viewport_height = 2.0 * h;
//    float viewport_width = aspect_ratio * viewport_height;
//
//    vec3 w = normalize(eye - lookat);
//    vec3 u = normalize(cross(up, w));
//    vec3 v = cross(w, u);
//
//    horizontal = viewport_width * u;
//    vertical = viewport_height * v;
//    lower_left_corner = eye - horizontal / 2.0f - vertical / 2.0f - w;
//}

//Ray Camera::get_ray(float s, float t)
//{
//
//    Ray ray;
//    ray.startPoint = eye;
//    ray.direction = lower_left_corner + s * horizontal + t * vertical - eye;
//    ray.direction = normalize(ray.direction);
//
//    return ray;
//}
//
//void Camera::Print(){
//    printf("Camera:\n");
//    printf("fovy: %f ", fovy);
//    printf("eye: (%f, %f, %f) ", eye.x, eye.y, eye.z);
//    printf("lookat: (%f, %f, %f) ", lookat.x, lookat.y, lookat.z);
//    printf("up: (%f, %f, %f) \n", up.x, up.y, up.z);
//}



//void Material::readinMap() {
//    std::cout << "Reading file " << map_Kd << std::endl;
//
//    img = cv::imread(map_Kd); //read in file
//
//    if (img.empty()) {
//        std::cout << "Cannot read file: " << map_Kd << std::endl;
//    }
//    map_height = img.rows, map_width = img.cols;
//}

//double Triangle::calAera()
//{
//    double a = length(p[1] - p[0]), b = length(p[2] - p[0]), c = length(p[2] - p[1]); // length of a,b,c
//    double cos_c = (a * a + b * b - c * c) / (2 * a * b);                             // 余弦定理
//    double sin_c = sqrt(1 - pow(cos_c, 2));
//    double aera = a * b * sin_c / 2;
//    return aera;
//}
//
//vec3 Triangle::findBaryCor(vec3 hitp)
//{
//    Eigen::Matrix<double, 4, 3> A;
//    Eigen::Matrix<double, 4, 1> B;
//    Eigen::MatrixXd res;
//
//    A << p[0].x, p[1].x, p[2].x,
//        p[0].y, p[1].y, p[2].y,
//        p[0].z, p[1].z, p[2].z,
//        1, 1, 1;
//    B << hitp.x, hitp.y, hitp.z, 1;
//
//    res = A.colPivHouseholderQr().solve(B);
//    vec3 ret;
//    ret.x = res(0, 0), ret.y = res(1, 0), ret.z = res(2, 0);
//
//    return ret;
//}