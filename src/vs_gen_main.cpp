#define _CRT_SECURE_NO_WARNINGS

#if defined(_WIN32) | defined(__WIN32__) | defined(__WIN32) | defined(_WIN64) | defined(__WIN64)
#include <windows.h>

#else

#endif

// C/C++ includes
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <type_traits>
#include <list>
#include <set>

// OpenCV includes
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

// custom includes
#include <cv_blur_process.h>
#include <cv_random_image_gen.h>
#include <cv_create_gaussian_kernel.h>
#include <cv_dft_conv.h>
#include <blur_params.h>
#include <num2string.h>
#include <file_ops.h>

//template<typename T>
//bool comp_pair(T fisrt, T second)
//{
//    return first < second;
//}


// ----------------------------------------------------------------------------------------
bool compare(std::pair<uint8_t, uint8_t> p1, std::pair<uint8_t, uint8_t> p2)
{
    return max(p1.first, p1.second) < max(p2.first, p2.second);
}


// ----------------------------------------------------------------------------------------
int main(int argc, char** argv)
{
    uint32_t idx = 0, jdx = 0;
    uint32_t img_h = 512;
    uint32_t img_w = 512;
    cv::Size img_size(img_h, img_w);

    cv::RNG rng(time(NULL));

    // timing variables
    typedef std::chrono::duration<double> d_sec;
    auto start_time = chrono::system_clock::now();
    auto stop_time = chrono::system_clock::now();
    auto elapsed_time = chrono::duration_cast<d_sec>(stop_time - start_time);
    std::string platform;

    cv::Mat img_f1, img_f2;
    cv::Mat kernel;
    cv::Mat output_img, mask;
    cv::Mat f1_layer, f2_layer;
    cv::Mat random_img;
    std::vector<uint16_t> dm_values;
    std::vector<uint16_t> dm_indexes;
    int32_t min_N, max_N;
    uint32_t N;
    double scale = 0.1;
    uint32_t BN = 1000;

    std::string scenario_name;
    std::pair<uint8_t, double> bg_dm;
    std::pair<uint8_t, double> fg_dm;
    std::vector<std::pair<uint8_t, uint8_t>> bg_br_table;
    std::vector<std::pair<uint8_t, uint8_t>> fg_br_table;
    double prob_bg = 0.31;    // set the probablility of selecting the background depthmap value
    double prob_fg = 0.35;    // set the probability of selecting the foreground depthmap value
    double bg_x = 0, fg_x = 0;

    std::vector<cv::Mat> blur_kernels;
    std::vector<cv::Mat> fft_blur_kernels;

    std::vector<uint8_t> depthmap_values;
    std::vector<double> sigma_table;
    std::vector<uint8_t> br1_table, tmp_br1_table;
    std::vector<uint8_t> br2_table, tmp_br2_table;
    
    // uint8_t aperture;
    // uint16_t slope;
    // uint16_t intercept;
    // uint32_t wavelength_min;
    // uint32_t wavelength_max;
    // double refractive_index_min;
    // double refractive_index_max;
    
    int32_t max_dm_num;
    uint32_t num_objects;


}   // end of main

