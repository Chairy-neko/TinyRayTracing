#include "svpng.inc"
#include "scene.h"
#include "bvh.h"
#include "pathtracing.h"

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <omp.h> // openmp accelerate
#include <time.h>
using namespace std;

int SAMPLE = 256;
string basedir = "";
string mtl_path = "";
string xml_path = "";
string obj_path = "";

void imshow(double *SRC, string index, int img_width, int img_height)
{

    unsigned char *image = new unsigned char[img_width * img_height * 3]; // img buffer
    unsigned char *p = image;
    double *S = SRC; //

    string name = basedir + "/image" + index + ".png";
    FILE *fp;
    fopen_s(&fp, name.c_str(), "wb");

    for (int i = 0; i < img_height; i++)
    {
        for (int j = 0; j < img_width; j++)
        {
            *p++ = (unsigned char)clamp(pow(*S++, 1.0f / 2.2f) * 255, 0.0, 255.0); // R
            *p++ = (unsigned char)clamp(pow(*S++, 1.0f / 2.2f) * 255, 0.0, 255.0); // G
            *p++ = (unsigned char)clamp(pow(*S++, 1.0f / 2.2f) * 255, 0.0, 255.0); // B
        }
    }

    svpng(fp, img_width, img_height, image, 0);
    cout << "\nIamge output to " << name << endl;
}

int main()
{
    printf("Please input base directory of the scene:\n");
    cin >> basedir;
    printf("Please input .mtl file path of the scene:\n");
    cin >> mtl_path;
    printf("Please input .xml file path of the scene:\n");
    cin >> xml_path;
    printf("Please input .obj file path of the scene:\n");
    cin >> obj_path;
    printf("Please input SPP:\n");
    scanf("%d", &SAMPLE);

    static std::default_random_engine e(time(NULL));
    static std::uniform_real_distribution<double> u1(0, 1);

    clock_t start, end, current;
    start = clock();
    Scene scene;
    set<int> remains;
    int tmp = 0;

    // load file, the order cannot be changed
    scene.readxml(xml_path);
    scene.readobj(obj_path);
    scene.readmtl(mtl_path, basedir);

    printf("image info:\nwidth: %d height: %d\n", scene.img_width, scene.img_height);
    scene.camera.Print();

    double *image = new double[scene.img_width * scene.img_height * 3];
    memset(image, 0.0, sizeof(double) * scene.img_width * scene.img_height * 3);
    BVHNode *root = buildBVH(scene.triangles, 0, scene.triangles.size() - 1, 8);
    printf("Build BVH down.\n");

    omp_set_num_threads(50); // parallel
#pragma omp parallel for
    for (int k = 0; k < SAMPLE; k++)
    {
        double *p = image;
        for (int i = 0; i < scene.img_height; i++)
        {
            for (int j = 0; j < scene.img_width; j++)
            {
                double x = double(j) / double(scene.img_width - 1.0);
                double y = double(scene.img_height - i) / double(scene.img_height - 1.0);

                // MSAA
                x += (u1(e) - 0.5f) / double(scene.img_width);
                y += (u1(e) - 0.5f) / double(scene.img_height);

                Ray ray = scene.camera.getRay(x, y);

                HitRecord rec = traverseBVH(ray, scene.triangles, root);
                vec3 color = vec3(0);
                if (rec.is_hit)
                {
                    color = shade(rec, -ray.direction, scene, root) / (float)SAMPLE;
                }
                *p += color.x;
                p++; // R
                *p += color.y;
                p++; // G
                *p += color.z;
                p++; // B
            }
            if (i % 100 == 0)
                printf("%d-%d ", k, i);
        }
    }
    imshow(image, to_string(SAMPLE), scene.img_width, scene.img_height);
    std::cerr << "\nDone.\n";
    end = clock();
    cout << (double)(end - start) / CLOCKS_PER_SEC << endl;

    return 0;
}