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
//#include <blur_params.h>
#include <num2string.h>
#include <file_ops.h>

#include <vs_gen.h>

//template<typename T>
//bool comp_pair(T fisrt, T second)
//{
//    return first < second;
//}


//-----------------------------------------------------------------------------
// global library internal state variables
vs_gen vs;


//-----------------------------------------------------------------------------
void init_from_file(std::string filename)
{


}

//-----------------------------------------------------------------------------
void init(/*params*/)
{

}


//-----------------------------------------------------------------------------
void generate_scene(unsigned int img_w, unsigned int img_h, unsigned char* img_ptr_f1, unsigned char* img_ptr_f2, unsigned short* dm_ptr)
{
    uint32_t idx, N;
    int32_t min_N, max_N;
    std::vector<uint8_t> tmp_br1_table, tmp_br2_table;
    std::vector<uint16_t> dm_values;
    cv::Mat random_img, output_img, mask;
    cv::Mat f1_layer, f2_layer;
    double scale2;

    // assigned the input pointers to cv:Mat containers
    cv::Mat img_f1 = cv::Mat(img_h, img_w, CV_8UC3, img_ptr_f1);
    cv::Mat img_f2 = cv::Mat(img_h, img_w, CV_8UC3, img_ptr_f2);
    cv::Mat dm = cv::Mat(img_h, img_w, CV_16UC1, dm_ptr);
    cv::Size img_size(img_h, img_w);

    // generate random dm_values that include the foreground and background values
    int32_t tmp_dm_num = vs.max_dm_vals_per_image;

    // get the probablility that the background depthmap value will be used
    double bg_x = vs.rng.uniform(0.0, 1.0);

    // get the probability that the foreground depthmap value will be used
    double fg_x = vs.rng.uniform(0.0, 1.0);

    if (bg_x < vs.prob_bg)
        tmp_dm_num--;

    if (fg_x < vs.prob_fg)
        tmp_dm_num--;

    //generate_depthmap_set(min_dm_value, max_dm_value, tmp_dm_num, depthmap_values, dm_values, rng);
    generate_depthmap_index_set(vs.min_dm_value, vs.max_dm_value, tmp_dm_num, vs.depthmap_values, vs.dm_indexes, vs.rng);

    // check the background probability and fill in the tables
    if (bg_x < vs.prob_bg)
    {
        uint16_t dm = vs.rng.uniform(0, vs.bg_br_table.size());
        tmp_br1_table.push_back(vs.bg_br_table[dm].first);
        tmp_br2_table.push_back(vs.bg_br_table[dm].second);

        dm_values.push_back(vs.bg_dm.first);
    }

    // fill in the tables for the region of interest depthmap values
    for (idx = 0; idx < vs.dm_indexes.size(); ++idx)
    {
        tmp_br1_table.push_back(vs.br1_table[vs.dm_indexes[idx]]);
        tmp_br2_table.push_back(vs.br2_table[vs.dm_indexes[idx]]);
        dm_values.push_back(vs.depthmap_values[vs.dm_indexes[idx]]);
    }

    // check the foreground probability and fill in the tables
    if (fg_x < vs.prob_fg)
    {
        uint16_t dm = vs.rng.uniform(0, vs.fg_br_table.size());
        tmp_br1_table.push_back(vs.fg_br_table[dm].first);
        tmp_br2_table.push_back(vs.fg_br_table[dm].second);

        dm_values.push_back(vs.fg_dm.first);
    }

    //N = (uint32_t)(img_h * img_w * 0.001);
    switch (vs.dataset_type)
    {
    case 0:
        generate_random_image(img_f1, vs.rng, img_h, img_w, vs.BN, vs.scale);
        break;

    //case 1:
    //    img_f1 = cv::Mat(img_h, img_w, CV_8UC3);
    //    sn_scale = sn_slope * dm_values[0] + sn_int;
    //    create_color_map(img_h, img_w, sn_scale, octaves, persistence, wood.data(), img_f1.data);
    //    break;

    //case 2:
    //    //img_f1 = cv::Mat(img_h, img_w, CV_8UC3, cv::Scalar::all(0,0,0));
    //    cv::hconcat(cv::Mat(img_h, img_w >> 1, CV_8UC3, cv::Scalar::all(0)), cv::Mat(img_h, img_w - (img_w >> 1), CV_8UC3, cv::Scalar::all(255)), img_f1);
    //    break;

    //    // black and white checkerboard
    //case 3:
    //    generate_checkerboard(48, 48, img_w, img_h, img_f1);
    //    break;
    }

    // clone the images
    img_f2 = img_f1.clone();

    //cv::Mat img_f1_t = img_f1.clone();
    //cv::Mat dst;
    //dft_conv_rgb(img_f1_t, fft_blur_kernels[tmp_br1_table[0]], dst);

    // create gaussian kernel and blur imgs
    cv::filter2D(img_f1, img_f1, -1, vs.blur_kernels[tmp_br1_table[0]], cv::Point(-1, -1), 0.0, cv::BorderTypes::BORDER_REPLICATE);
    cv::filter2D(img_f2, img_f2, -1, vs.blur_kernels[tmp_br2_table[0]], cv::Point(-1, -1), 0.0, cv::BorderTypes::BORDER_REPLICATE);

    // create the initial depth map
    cv::Mat depth_map(img_h, img_w, CV_8UC1, cv::Scalar::all(dm_values[0]));

    // blur imgs using dm_values and random masks
    for (idx = 1; idx < dm_values.size(); ++idx)
    {
        f1_layer = img_f1.clone();
        f2_layer = img_f2.clone();

        // This part help shape the final distribution of depthmap values
        //min_N = (int32_t)ceil((num_objects) / (1 + exp(-0.35 * dm_values[idx] + (0.035 * num_objects))) + 2);
        //min_N = (int32_t)(num_objects / (double)(1.0 + exp(-0.1 * (depthmap_values[dm_indexes[idx]] - (max_dm_value-min_dm_value)/2.0))) + 5);
        //min_N = (int32_t)ceil(((max_dm_value) / (double)(1.0 + exp(-0.35 * depthmap_values[dm_values[idx]] + (0.175 * max_dm_value)))) + 3);
        min_N = (int32_t)ceil(((vs.max_dm_value) / (double)(1.0 + exp(-0.365 * dm_values[idx] + (0.175 * vs.max_dm_value)))) + 3);
        max_N = (int32_t)ceil(2.0 * min_N);  // 2.0

        N = vs.rng.uniform(min_N, max_N + 1);

        // define the scale factor
        scale2 = 60.0 / (double)img_size.width;

        // 
        switch (vs.dataset_type)
        {
        case 0:
            generate_random_image(random_img, vs.rng, img_h, img_w, vs.BN, scale2);
            break;

        //case 1:
        //    sn_scale = sn_slope * dm_values[idx] + sn_int;
        //    random_img = cv::Mat(img_h, img_w, CV_8UC3);
        //    create_color_map(img_h, img_w, sn_scale, octaves, persistence, wood.data(), random_img.data);
        //    break;

        //case 2:
        //    cv::hconcat(cv::Mat(img_h, img_w >> 1, CV_8UC3, cv::Scalar::all(0)), cv::Mat(img_h, img_w - (img_w >> 1), CV_8UC3, cv::Scalar::all(255)), random_img);
        //    break;
        }

        // generate random overlay
        generate_random_overlay(random_img, vs.rng, output_img, mask, N, vs.scale);

        overlay_image(f1_layer, output_img, mask);
        overlay_image(f2_layer, output_img, mask);

        // overlay depthmap
        overlay_depthmap(depth_map, mask, dm_values[idx]);

        // blur f1
        blur_layer(f1_layer, img_f1, mask, vs.blur_kernels[tmp_br1_table[idx]], vs.rng);

        // blur f2
        blur_layer(f2_layer, img_f2, mask, vs.blur_kernels[tmp_br2_table[idx]], vs.rng);
    }

}   // end of generate_scene




////-----------------------------------------------------------------------------
//int main(int argc, char** argv)
//{
//    uint32_t idx = 0, jdx = 0;
//    uint32_t img_h = 512;
//    uint32_t img_w = 512;
//    cv::Size img_size(img_h, img_w);
//
//    cv::RNG rng(time(NULL));
//
//    // timing variables
//    typedef std::chrono::duration<double> d_sec;
//    auto start_time = chrono::system_clock::now();
//    auto stop_time = chrono::system_clock::now();
//    auto elapsed_time = chrono::duration_cast<d_sec>(stop_time - start_time);
//    std::string platform;
//
//    cv::Mat img_f1, img_f2;
//    cv::Mat kernel;
//    cv::Mat output_img, mask;
//    cv::Mat f1_layer, f2_layer;
//    cv::Mat random_img;
//    std::vector<uint16_t> dm_values;
//    std::vector<uint16_t> dm_indexes;
//    int32_t min_N, max_N;
//    uint32_t N;
//    double scale = 0.1;
//    uint32_t BN = 1000;
//
//    std::string scenario_name;
//    std::pair<uint8_t, double> bg_dm;
//    std::pair<uint8_t, double> fg_dm;
//    std::vector<std::pair<uint8_t, uint8_t>> bg_br_table;
//    std::vector<std::pair<uint8_t, uint8_t>> fg_br_table;
//    double prob_bg = 0.31;    // set the probablility of selecting the background depthmap value
//    double prob_fg = 0.35;    // set the probability of selecting the foreground depthmap value
//    double bg_x = 0, fg_x = 0;
//
//    std::vector<cv::Mat> blur_kernels;
//    std::vector<cv::Mat> fft_blur_kernels;
//
//    std::vector<uint8_t> depthmap_values;
//    std::vector<double> sigma_table;
//    std::vector<uint8_t> br1_table, tmp_br1_table;
//    std::vector<uint8_t> br2_table, tmp_br2_table;
//    
//    // uint8_t aperture;
//    // uint16_t slope;
//    // uint16_t intercept;
//    // uint32_t wavelength_min;
//    // uint32_t wavelength_max;
//    // double refractive_index_min;
//    // double refractive_index_max;
//    
//    int32_t max_dm_num;
//    uint32_t num_objects;
//
//
//}   // end of main

