#pragma once
#include "../dataStructure/BVH.hpp"
#include "../dataStructure/Ray.hpp"
#include <queue>
#include <cmath>
#include <omp.h>
#include <random>
#include <glm/glm.hpp>
#include "svpng.inc"

static omp_lock_t lock;
bool rayIntersectScene(Ray& ray, BVH& B, Intersection& p) {
	Intersection i;
	std::queue<BVNode> bvq;
	if (ray.intersectBoundingbox(B.bvh[1].box)) {
		bvq.push(B.bvh[1]);
	}
	while (!bvq.empty()) {
		BVNode node = bvq.front();
		bvq.pop();
		if (node.leftchild == -1 && node.rightchild == -1) { //Ҷ�ڵ�
			Vertex temp;
			if (ray.intersectFace(*node.object, temp)) {
				double t = (temp - ray.start).norm() / ray.direction.norm();
				if ((temp - ray.start).dot(ray.direction) < 0) {
					t = -t;
				}
				if (t > 0 && (!i.valid || (i.valid && t < i.t))) {
					i.valid = true;
					i.t = t;
					i.ray = ray;
					i.f = *node.object;
					i.p = temp;
				}
			}

		}
		else {
			int leftchild = node.leftchild;
			if (leftchild > 0 && ray.intersectBoundingbox(B.bvh[leftchild].box)) {
				bvq.push(B.bvh[leftchild]);
			}
			int rightchild = node.rightchild;
			if (rightchild > 0 && ray.intersectBoundingbox(B.bvh[rightchild].box)) {
				bvq.push(B.bvh[rightchild]);
			}
		}
	}
	if (i.valid) {
		Eigen::Vector3d bc = util.barycentric_coordinate(i.f, i.p);
		i.pn = (i.f.vn1 * bc[0] + i.f.vn2 * bc[1] + i.f.vn3 * bc[2]).normalized();
		p = i;
		return true;
	}
	else {
		return false;
	}

}


bool russian_Roulette(double pr)
{
	static std::default_random_engine e;
	static std::uniform_real_distribution<double>u1(0, 1);
	double rnd = u1(e);
	if (rnd < pr) return true;
	else return false;
}

bool Refract(Eigen::Vector3d dir, Eigen::Vector3d& normal, double eta, Vertex& refract_dir) {
	Vertex& i = dir;
	Vertex& n = normal;
	float cosi = i.dot(n);
	float cost2 = 1.0f - eta * eta * (1.0f - cosi * cosi);
	if (cost2 >= 0.0f) {
		refract_dir = i * eta - n * (eta * cosi + std::sqrt(cost2));
		return true;
	}
	else {
		return false;
	}
}

Vertex BRDFImportanceSampling(Vertex& direction, rayType type, Face& f, double Ns) {

	static std::default_random_engine e(time(NULL));
	static std::uniform_real_distribution<double>u1(0, 1);

	double phi = u1(e) * 2 * PI;
	double theta;
	if (type == DIFFUSE) {
		theta = asin(sqrt(u1(e)));
	}
	else if (type == SPECULAR) {
		theta = acos(pow(u1(e), (double)1 / (Ns + 1)));
	}

	Vertex sample(sin(theta) * cos(phi), cos(theta), sin(theta) * sin(phi));
	Vertex front;
	if (fabs(direction[0]) > fabs(direction[1])) {
		front = Vertex(direction[2], 0, -direction[0]).normalized();
	}
	else {
		front = Vertex(0, -direction[2], direction[1]).normalized();
	}
	Vertex right = direction.cross(front);
	Vertex ret = (right * sample[0]) + (direction * sample[1]) + (front * sample[2]);
	ret.normalize();
	return ret;
}

Ray nextRay(Intersection& p, Eigen::Vector3d dir, Scene& scene) {
	static std::default_random_engine e(time(NULL));
	static std::uniform_real_distribution<double>u2(0, 1);

	Material* material = scene.material_map[p.f.material];
	Eigen::Vector3d direction;
	//͸������
	if (material->Ni > 1) {
		double n1, n2;
		double cos_in = -dir.dot(p.pn);
		Vertex normal;
		if (cos_in > 0) {
			normal = -p.pn;
			n1 = material->Ni; //�������
			n2 = 1.0;
		}
		else {
			normal = p.pn;
			n1 = 1.0;
			n2 = material->Ni; //�������
		}

		double rf0 = pow((n1 - n2) / (n1 + n2), 2);
		double fresnel = rf0 + (1.0f - rf0) * pow(1.0f - std::abs(cos_in), 5);
		if (fresnel < u2(e)) {
			if (Refract(-dir, normal, n1 / n2, direction)) {  //������ҵ����䷽��
				return Ray(p.p, direction, TRANSMISSION);  //����һ�������䷽��Ĺ���
			}
			else {
				Eigen::Vector3d reflect;
				Eigen::Vector3d incoming = -dir;   //�����ؾ��淴�䷽������һ���������
				reflect = incoming - normal * (incoming.dot(normal)) * 2;
				return Ray(p.p, reflect, SPECULAR);
			}
		}
	}
	//��͸������
	Eigen::Vector3d kd = material->kd, ks = material->ks;
	double kd_av = (kd[0] + kd[1] + kd[2]) / 3.0, ks_av = (ks[0] + ks[1] + ks[2]) / 3.0;
	rayType ray_type;
	double prob = u2(e);
	if (prob < kd_av) {
		direction = BRDFImportanceSampling(p.pn, DIFFUSE, p.f, material->Ns);
		ray_type = DIFFUSE;
	}
	else if (prob < kd_av + ks_av) {
		Eigen::Vector3d incoming = -dir;
		Eigen::Vector3d reflect = incoming - p.pn * (incoming.dot(p.pn)) * 2;  //���淴�䣬���ķ���Ϊ���淴�䷽��
		direction = BRDFImportanceSampling(reflect, SPECULAR, p.f, material->Ns);
		ray_type = SPECULAR;
	}
	else {
		ray_type = INVALID;
	}

	Ray r(p.p + (direction * 1e-3), direction, ray_type);
	return r;
}

Eigen::Vector3d shade(Intersection& p, Eigen::Vector3d dir, Scene& scene, BVH& B)
{
	if (scene.light_map.find(p.f.material) != scene.light_map.end()) {  //������ߴ��۾��������������غ����̴��˹�Դ��
		return scene.light_map[p.f.material]->radiance;  //ֱ�ӷ��ع�Դ�ķ�����ɫ
	}
	//����������ϣ����ȼ���򵽵ĵ�p���Ĳ���ϵ��kd
	Material* material = scene.material_map[p.f.material];
	Eigen::Vector3d kd, ks;
	if (material->mapping_flag == true) { //����˴�������
		Eigen::Vector3d bc = util.barycentric_coordinate(p.f, p.p);
		Vt vt = p.f.vt1 * bc[0] + p.f.vt2 * bc[1] + p.f.vt3 * bc[2];
		vt[0] = vt[0] - floor(vt[0]); vt[1] = vt[1] - floor(vt[1]);
		int vt0 = vt[0] * material->map_height, vt1 = vt[1] * material->map_width;
		cv::Vec3b vec_3 = material->texture_map.at<cv::Vec3b>(vt0, vt1); //��ȡ����ֵ
		kd[0] = (double)vec_3[2] / 255, kd[1] = (double)vec_3[1] / 255, kd[2] = (double)vec_3[0] / 255;  //opencv����ͼƬ����BGR��˳���Ÿ���ɫͨ����ֵ����������ǵ�˳�������RGB
	}
	else {
		kd = material->kd;
	}
	ks = material->ks;
	double Ns = material->Ns;

	Eigen::Vector3d L_dir(0, 0, 0), L_indir(0, 0, 0);
	static std::default_random_engine e(time(NULL));
	//ֱ�ӹ���
	for (int i = 0; i < scene.l.size(); i++) {
		Light l = scene.l[i];
		std::vector<double> accu_area = scene.material_map[l.name]->accum_area;  //����Դ�������������
		double total_area = accu_area[accu_area.size() - 1];
		static std::uniform_real_distribution<double>u1(0, total_area);
		double rnd = u1(e);
		int sample = 0;
		while (accu_area[sample] <= rnd && sample < accu_area.size() - 1) {
			sample++;  //
		}
		Face sample_face = scene.material_map[l.name]->f[sample];  //����������������Ƭ
		static std::uniform_real_distribution<double>u2(0, 1);   //������һ���������꣬�������������ɸ����������ϵ�һ����
		double rnd1 = u2(e), rnd2 = u2(e), rnd3 = u2(e);
		double p1 = rnd1 / (rnd1 + rnd2 + rnd3), p2 = rnd2 / (rnd1 + rnd2 + rnd3), p3 = rnd3 / (rnd1 + rnd2 + rnd3);
		Vertex xl = sample_face.v1 * p1 + sample_face.v2 * p2 + sample_face.v3 * p3;  //�������Ĺ�Դ��
		Vn vn = (sample_face.vn1 * p1 + sample_face.vn2 * p2 + sample_face.vn3 * p3).normalized();  //�õ㴦��Դ��ķ���

		Eigen::Vector3d direction = (xl - p.p).normalized(); //�Ӹõ㵽��Դ����������߷���
		double visibility = 1;
		Ray rl;
		rl.start = p.p + (direction * 1e-3);//add micro turbulance
		rl.direction = direction;
		Intersection inter;
		rayIntersectScene(rl, B, inter);
		if (inter.f.material != sample_face.material) {
			visibility = 0;
		}

		if (direction.dot(p.pn) > 0) {
			double pdf_light = double(1) / total_area;
			double cos_theta = abs(direction.dot(vn));
			double cos_theta_hat = abs(direction.dot(p.pn) / p.pn.norm());
			Eigen::Vector3d intensity = l.radiance * cos_theta * cos_theta_hat / pow((xl - p.p).norm(), 2) / pdf_light * visibility;
			
			Eigen::Vector3d reflect_direction = 2 * (direction.dot(p.pn)) * p.pn - direction;
			double cos_alpha = fmax(dir.dot(reflect_direction), 0.0);
			L_dir[0] += intensity[0] * (kd[0] / PI  + ks[0] * (Ns + 2) * pow(cos_alpha, Ns) / (2 * PI));
			L_dir[1] += intensity[1] * (kd[1] / PI  + ks[1] * (Ns + 2) * pow(cos_alpha, Ns) / (2 * PI));
			L_dir[2] += intensity[2] * (kd[2] / PI  + ks[2] * (Ns + 2) * pow(cos_alpha, Ns) / (2 * PI));
			
		}
	}
	//��ӹ���
	double P_RR = 0.5; //����˹���̶ĵĸ���
	if (russian_Roulette(P_RR)) {
		Ray r = nextRay(p, dir, scene); //�������ɹ���dir��ǰһ�����ߣ� ��BRDF����Ҫ�Բ����õ�
		Intersection ret;
		if (r.ray_type != INVALID && rayIntersectScene(r, B, ret)) {
			Vertex intensity = shade(ret, -1 * r.direction, scene, B) / P_RR;
			double cos_theta_hat = abs(r.direction.dot(p.pn) / p.pn.norm());
			if (r.ray_type == DIFFUSE) {//�����޹������壬��Ϊ������
				if (scene.light_map.find(ret.f.material) == scene.light_map.end()) {
					L_indir[0] +=  intensity[0] * kd[0];
					L_indir[1] +=  intensity[1] * kd[1];
					L_indir[2] +=  intensity[2] * kd[2];
				}
			}
			else if (r.ray_type == SPECULAR) {
				if (scene.light_map.find(ret.f.material) == scene.light_map.end()) {
					Eigen::Vector3d reflect_direction = 2 * (r.direction.dot(p.pn)) * p.pn - r.direction;
					double cos_alpha = dir.dot(reflect_direction);
					L_indir[0] += intensity[0] * ks[0];
					L_indir[1] += intensity[1] * ks[1];
					L_indir[2] += intensity[2] * ks[2];
				}
			}
			else {
				L_indir = L_indir + intensity;
			}
		}
	}
	return L_dir + L_indir;
}

void generateImg(Scene& scene, BVH& B, Image& img, int N_ray_per_pixel) {

	Eigen::Vector3d dir = scene.camera.look_at - scene.camera.eye;
	double l = dir.norm();
	double screen_height = 2 * l * tan(scene.camera.fovy * PI / 360);
	double screen_width = screen_height / scene.camera.height * scene.camera.width;

	Vertex screen_center = scene.camera.look_at; //��Ļ���ĵ����������
	double pdx = screen_width / scene.camera.width, pdy = screen_height / scene.camera.height; //ÿ�������������еĴ�С

	scene.camera.up.normalize();
	Eigen::Vector3d screen_y_dir = scene.camera.up; //ȡ������ϵķ���Ϊy�������򣬿���ķ���Ϊz��������������ҵķ���Ϊx����������x���� = z���� x y����
	Eigen::Vector3d screen_x_dir = dir.normalized().cross(screen_y_dir).normalized();

	Eigen::Vector3d screen_pdy = screen_y_dir * pdy; //��Ļ��y�������ƶ�һ�����ص�λ��ʸ��������������ϵ�µ�λ��
	Eigen::Vector3d screen_pdx = screen_x_dir * pdx; //��Ļ��x�������ƶ�һ�����ص�λ��ʸ��������������ϵ�µ�λ��

	//����Ļ�����Ͻǿ�ʼö������
	Eigen::Vector3d left_top_corner = screen_center - 0.5 * screen_width * screen_x_dir + 0.5 * screen_height * screen_y_dir;
	Eigen::Vector3d screen_pos;
	omp_init_lock(&lock);

	//for (int i = 0; i < scene.camera.height; i += 2) { //���ϵ���
	//	for (int j = 0; j < scene.camera.width; j += 2) { //
	for (int i = 0; i < scene.camera.height; i++) { //���ϵ���
		for (int j = 0; j < scene.camera.width; j++) { //������
			screen_pos = left_top_corner - i * screen_pdy + j * screen_pdx;
#pragma omp parallel for
			for (int k = 0; k < N_ray_per_pixel; k++) {
				Ray ray;
				ray.start = scene.camera.eye;
				ray.direction = (screen_pos - ray.start).normalized();
				//�����������볡���Ľ���
				Intersection p;
				if (rayIntersectScene(ray, B, p)) {
					Eigen::Vector3d radiance = shade(p, -ray.direction, scene, B);
					omp_set_lock(&lock);
					img.img[img.getIndex(j, i, 0)] += radiance[0] / N_ray_per_pixel;
					img.img[img.getIndex(j, i, 1)] += radiance[1] / N_ray_per_pixel;
					img.img[img.getIndex(j, i, 2)] += radiance[2] / N_ray_per_pixel;
					/*img.img[img.getIndex(j + 1, i, 0)] += radiance[0] / N_ray_per_pixel;
					img.img[img.getIndex(j + 1, i, 0)] += radiance[1] / N_ray_per_pixel;
					img.img[img.getIndex(j + 1, i, 0)] += radiance[2] / N_ray_per_pixel;
					img.img[img.getIndex(j, i + 1, 0)] += radiance[0] / N_ray_per_pixel;
					img.img[img.getIndex(j, i + 1, 1)] += radiance[1] / N_ray_per_pixel;
					img.img[img.getIndex(j, i + 1, 2)] += radiance[2] / N_ray_per_pixel;
					img.img[img.getIndex(j + 1, i + 1, 0)] += radiance[0] / N_ray_per_pixel;
					img.img[img.getIndex(j + 1, i + 1, 1)] += radiance[1] / N_ray_per_pixel;
					img.img[img.getIndex(j + 1, i + 1, 2)] += radiance[2] / N_ray_per_pixel;*/
					omp_unset_lock(&lock);
				}
			}
		}
		if (i % 16 == 0) {
			std::cout << "The " << i << "-th row finish rendering.\n";
		}
	}
	omp_destroy_lock(&lock);
}

void saveImg(double* SRC, int height, int width, std::string filename, int N_ray_per_pixel) {

	unsigned char* image = new unsigned char[height * width * 3];// ͼ��buffer
	unsigned char* p = image;
	double* S = SRC;    // Դ����

	FILE* fp;
	char* buffer = new char[10];
	_itoa(N_ray_per_pixel, buffer, 10);
	fopen_s(&fp, (filename + "-SPP" + buffer + ".png").c_str(), "wb");

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			*p++ = (unsigned char)glm::clamp(pow((*S++), 0.45454545) * 255, 0.0, 255.0);  // R ͨ��
			*p++ = (unsigned char)glm::clamp(pow((*S++), 0.45454545) * 255, 0.0, 255.0);  // G ͨ��
			*p++ = (unsigned char)glm::clamp(pow((*S++), 0.45454545) * 255, 0.0, 255.0);  // B ͨ��
		}
	}
	svpng(fp, width, height, image, 0);
}