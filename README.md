# Document

## 开发环境

- 操作系统：Windows 11 x64
- 处理器：12th Gen Intel(R) Core(TM) i9-12900H   2.50 GHz
- IDE：Microsoft Visual Studio 2019 x64
- 运行环境：Release X64
- 依赖库：glm（向量运算），Eigen（高性能数学运算），opencv（texture贴图），tinyxml2（xml文件读取）

## 使用说明

打开.exe文件，按照提示输入对应文件地址。生成结果会输出图像地址和运行时间（单位：s）。

![image-20230331013130045](D:\Data\codes\vs\RayTracingOnCPU\README.assets\image-20230331013130045.png)

## 数据结构

### Camera

相机类，存储相机信息。成员函数可通过相机信息设定机位，通过图片像素位置获取视线，以及进行相机信息的打印。

```c++
class Camera
{
public:
    Camera() {}
    ~Camera() {}

    void setCamera();// set camera 
    Ray getRay(float u, float v);// get ray at (u,v)
    void Print();// print info for debugging

    double fovy = 90;
    vec3 eye = vec3(278.0, 273.0, -800.0);
    vec3 lookat = vec3(278.0, 273.0, -799.0);
    vec3 up = vec3(0.0, 1.0, 0.0);
    double aspect_ratio = 1.0;

    vec3 lower_left_corner = vec3(0.0, 0.0, 0.0);
    vec3 horizontal = vec3(0.0, 0.0, 0.0);
    vec3 vertical = vec3(0.0, 0.0, 0.0);
};
```

### Material

Material类存储mtl读入的材质信息，`readinMap()`用于读入纹理地址并存储纹理数据

```c++
class Material
{
public:
    Material() {}
    ~Material() {}
    void readinMap(); //for texture

    vec3 Kd = vec3(0.0, 0.0, 0.0); // diffuse
    vec3 Ks = vec3(0.0, 0.0, 0.0); // specular
    vec3 Tr = vec3(0.0, 0.0, 0.0); // transmittance
    float Ns = 1;                  // shiness, the exponent of phong lobe
    float Ni = 1;                  // the Index of Refraction (IOR) of transparent object
    string map_Kd = "";            // map_Kd

    bool is_emissive = false;
    vec3 radiance = vec3(0, 0, 0);
    double area = 0.0;

    vector<Triangle> triangles; // for light sampling

    cv::Mat img;// for texture
    int map_height, map_width;
};
```

### Light

Light类是为了方便光源采样时遍历光源。

#### Light类

```c++
class Light
{
public:
    Light() {}
    Light(string m, vec3 r) : mtl_name(m), radiance(r) {}
    ~Light() {}

    string mtl_name = "";
    vec3 radiance = vec3(0.0, 0.0, 0.0);
};
```

### Ray

ray记录了一条射线的起点、终点和类型。

```c++
const int DIFFUSE      = 0;
const int SPECULAR     = 1;
const int TRANSMISSION = 2;
const int INVALID = 3;

class Ray
{
public:
    Ray() {}
    Ray(vec3 s, vec3 d) :startpoint(s), direction(d) {}
    Ray(vec3 s, vec3 d, int r) :startpoint(s), direction(d), ray_type(r) {}
    ~Ray() {}

    vec3 startpoint = vec3(0.0, 0.0, 0.0);
    vec3 direction = vec3(0.0, 0.0, 0.0);
    int ray_type = INVALID;
};
```

### Triangle

存储三角面片信息。Triangle类的成员函数可以计算三角形面积和取重心坐标。

```
class Triangle
{
public:
    double calAera();            // calulate area of the triangle
    vec3 findBaryCor(vec3 hitp); // fing the barycenter coordinate of the hitpoint

    Triangle() {}

    vec3 v[3];                         // vertex
    vec3 vn[3];                        // vertex normal
    vec2 vt[3];                        // vertex texture
    vec3 normal = vec3(0.0, 0.0, 0.0); // for caculating hitpoints and distances
    vec3 center = vec3(0.0, 0.0, 0.0); // for sorting
    double area = 0.0;                 // for light sampling
    std::string mtl_name = "";
    bool is_emissive = false; // choose Emissive triangle when they are overlapping
    // int id = 0; // for debug
};
```

### Scene

Scene类存储场景内的Material、Triangle、Light、Camera信息，成员函数对obj/xml/mtl文件进行读取。

```c++
class Scene
{
public:
    Scene(/* args */) {}
    ~Scene() {}

    void readxml(string xml_path);                  // get camera & lights
    void readmtl(string mtl_path, string base_dir); // get materials
    void readobj(string obj_path);                  // get triangles

    int img_width;
    int img_height;

    vector<Triangle> triangles; // obj models
    vector<Light> lights;
    unordered_map<string, Material> materials;
    Camera camera;
};
```

## 算法实现

#### 求交加速结构：基于SAH构建的BVH树

参考[2]基于SAH（Surface Area Heuristic，即表面积启发算法）建树。SAH 算法首先定义了查找一颗 BVH 树的代价为左子树查找时间和概率之积加右子树查找时间和概率之积，其中光线击中两个盒子的概率和盒子的表面积成正相关。

```c++
BVHNode *buildBVH(std::vector<Triangle> &triangles, int left_index, int right_index, int leaf_num)
{
    if (left_index > right_index)
        return 0;

    BVHNode *node = new BVHNode();
    node->AA = vec3(1145141919, 1145141919, 1145141919);
    node->BB = vec3(-1145141919, -1145141919, -1145141919);
    // caculate AABB
    for (int i = left_index; i <= right_index; i++)
    {
        // minimum point AA
        ...
        // maximum point BB
        ...
    }
    // if number of triangles < num_node, arrive at the leaf node
    if ((right_index - left_index + 1) <= leaf_num)
    {
        ...
    }
    // find the best split according to cost
    float Cost = INF;
    int Axis = 0;
    int Split = (left_index + right_index) / 2;
    for (int axis = 0; axis < 3; axis++)
    {
        // sort according to x/y/z axis
       	...
        // caculate cost
        float cost = INF;
        int split = left_index;
        for (int i = left_index; i <= right_index - 1; i++)
        {
            ...
            float leftCost = left_surface_area * (i - left_index + 1);
            float rightCost = right_surface_area * (right_index - i);
            float totalCost = leftCost + rightCost;
            if (totalCost < cost)
            {
                cost = totalCost;
                split = i;
            }
        }
        if (cost < Cost)
        {
            Cost = cost;
            Axis = axis;
            Split = split;
        }
    }
    // build BVH according to best Axis
    ...
    // build left node & right node
    node->left = buildBVH(triangles, left_index, Split, leaf_num);
    node->right = buildBVH(triangles, Split + 1, right_index, leaf_num);

    return node;
}
```

遍历BVH，击中三角形则返回该三角型，否则返回空指针。光线如果遇到叶子节点，和叶子节点的所有三角形求交。如果是中间节点，和左右子树的 AABB 盒子求交，如果击中盒子，就递归左右子树，取距离最近的求交结果，否则停止。

#### 蒙特卡洛路径追踪

代码逻辑参考games101[5]光线追踪部分伪代码。

<img src="D:\Data\codes\vs\RayTracingOnCPU\README.assets\image-20230331012028287.png" alt="image-20230331012028287" style="zoom:50%;" />

#### 直接光照

参考[4]中phong模型编写。

```c++
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
                ...
                bool visibility = true;                      // judge visibility
                ...
                HitRecord rec_sample = traverseBVH(ray_light, scene.triangles, root);
                if (rec_sample.triangle.mtl_name != light_triangle.mtl_name)
                {
                    visibility = false;
                }
                if (visibility && dot(wo, rec.pn) > 0)
                {
                    ...
                    vec3 intensity = radiance * cos_theta_p * cos_theta / length2(light_p - rec.hitpoint) / pdf_light;
                    ...
                    L_dir += intensity * (Kd / PI + m.Ks * (m.Ns + 2.0f) * (float)pow(cos_alpha, m.Ns) / (2.0f * PI));
                }
                break;
            }
        }
    }
```

#### 间接光照

参考[4]用多重重要性采样方法采样下一条光线，用俄罗斯轮盘赌控制递归次数（P_RR=0.8)。

```c++
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
```

## 测试结果

用openmp（50线程）在Release X64环境下测试结果如下。

### cornell-box

256spp, 72973.8s $\approx$ 20.3h (跑的时候忘了设置电脑不休眠了，时间可能有点问题)

<img src="D:\Data\codes\vs\RayTracingOnCPU\README.assets\image256.png" alt="image256" style="zoom:50%;" />

### staircase

256spp, 6088.5s $\approx$ 1.7h

<img src="D:\Data\codes\vs\RayTracingOnCPU\README.assets\image256-1680188181683-5.png" alt="image256" style="zoom:50%;" />

分析：天花板光和玻璃反射不对，整体场景偏亮。调了很久，分析这个场景和另外两个场景相比有更多光源，并且每种材质的光源不止有一处光源，所有点光源采样可能存在问题，导致漫反射和镜面反射的强度与方向与参考答案有较大差距。

### veach-mis

256spp, 268.1s $\approx$ 4.5min

<img src="D:\Data\codes\vs\RayTracingOnCPU\README.assets\image256-1680158811354-3.png" alt="image256" style="zoom:50%;" />

### own scene

10spp，主要是测试用的

<img src="D:\Data\codes\vs\RayTracingOnCPU\README.assets\image10.png" alt="image10" style="zoom:50%;" />

## 参考

[1] [_Ray Tracing in One Weekend_](https://raytracing.github.io/books/RayTracingInOneWeekend.html)

[2] [_光线追踪渲染实战（二）：BVH 加速遍历结构_](https://blog.csdn.net/weixin_44176696/article/details/118655688)

[3] [_tinyxml2读写XML文件的例程_](https://blog.csdn.net/zhawk/article/details/60880036)

[4] [_蒙特卡洛路径跟踪 C++实现_](https://blog.csdn.net/Listoree/article/details/124081645)

[4] Lafortune E P, Willems Y D. Using the modified phong reflectance model for physically based rendering[M]. Katholieke Universiteit Leuven. Departement Computerwetenschappen, 1994.

[5] [_GAMES101: 现代计算机图形学入门_](https://sites.cs.ucsb.edu/~lingqi/teaching/games101.html)
