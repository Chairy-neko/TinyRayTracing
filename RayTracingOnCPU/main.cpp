#include "svpng.inc"
#include "scene.h"
#include "bvh1.h"
#include "pathtracing.h"

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <omp.h> // openmp accelerate
#include <time.h>
using namespace std;

 const string basedir = "example-scenes-cg22/staircase";
 const string mtl_path = "example-scenes-cg22/staircase/staircase.mtl";
 const string xml_path = "example-scenes-cg22/staircase/staircase.xml";
 const string obj_path = "example-scenes-cg22/staircase/staircase.obj";
//const string basedir = "example-scenes-cg22/cornell-box";
//const string mtl_path = "example-scenes-cg22/cornell-box/cornell-box.mtl";
//const string xml_path = "example-scenes-cg22/cornell-box/cornell-box.xml";
//const string obj_path = "example-scenes-cg22/cornell-box/cornell-box.obj";
//const string basedir = "example-scenes-cg22/veach-mis";
//const string mtl_path = "example-scenes-cg22/veach-mis/veach-mis.mtl";
//const string xml_path = "example-scenes-cg22/veach-mis/veach-mis.xml";
//const string obj_path = "example-scenes-cg22/veach-mis/veach-mis.obj";
//const string basedir = "example-scenes-cg22/test";
//const string mtl_path = "example-scenes-cg22/test/back.mtl";
//const string xml_path = "example-scenes-cg22/test/back.xml";
//const string obj_path = "example-scenes-cg22/test/back.obj";

const int SAMPLE = 10;

void imshow(double *SRC, string index, int img_width, int img_height)
{

    unsigned char *image = new unsigned char[img_width * img_height * 3]; // img buffer
    // memset(image, 0.0, sizeof(double) * img_width * img_height * 3);
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
    BVHNode *root = buildBVHwithSAH(scene.triangles, 0, scene.triangles.size() - 1, 8);
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

                Ray ray = scene.camera.get_ray(x, y);

                HitResult res = hitBVH(ray, scene.triangles, root);
                vec3 color = vec3(0);
                if (res.isHit)
                {
                    color = shade(res, -ray.direction, scene, root) / (float)SAMPLE;
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
    printf("num_kd: %d, num_ks: %d", num_kd, num_ks);
    std::cerr << "\nDone.\n";
    end = clock();
    cout << (double)(end - start) / CLOCKS_PER_SEC << endl;

    return 0;
}