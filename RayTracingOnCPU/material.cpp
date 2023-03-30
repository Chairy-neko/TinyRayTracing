#include "material.h"

void Material::readinMap() {
    std::cout << "Reading map_Kd file " << map_Kd << std::endl;

    img = cv::imread(map_Kd);
    if (img.empty()) {
        std::cout << "Cannot read file: " << map_Kd << std::endl;
    }
    map_height = img.rows, map_width = img.cols;
}