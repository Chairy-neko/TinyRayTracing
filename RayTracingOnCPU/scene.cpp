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
    camera.setCamera();
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
        materials[mtl_name].is_emissive = true;
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
                        triangle.v[fi] = vertices[tmpindex];
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
            triangle.normal = normalize(cross(triangle.v[1] - triangle.v[0], triangle.v[2] - triangle.v[0]));
            triangle.center = (triangle.v[0] + triangle.v[1] + triangle.v[2]) / vec3(3.0, 3.0, 3.0);
            triangle.mtl_name = mtl_name;
            if(materials[mtl_name].is_emissive){
                triangle.is_emissive = true;
                double triangle_area = triangle.calAera();
                materials[mtl_name].area += triangle_area;// mtl total area
                triangle.area = materials[mtl_name].area;
                materials[mtl_name].triangles.push_back(triangle);
            }
            triangles.push_back(triangle);
        }
    }
    printf("num of vertices: %d\n", (int)vertices.size());
    printf("num of vn: %d\n", (int)vn.size());
    printf("num of vt: %d\n", (int)vt.size());
    printf("num of triangles: %d\n", (int)triangles.size());
}
