#include "material.h"

void Material::readinMap() {
    std::cout << "Reading file " << map_Kd << std::endl;

    img = cv::imread(map_Kd); //read in file

    if (img.empty()) {
        std::cout << "Cannot read file: " << map_Kd << std::endl;
    }
    map_height = img.rows, map_width = img.cols;
}