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
//#include <array>
#include <algorithm>
//#include <type_traits>
//#include <list>
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
//#include <num2string.h>
//#include <file_ops.h>

#include <vs_gen.h>
#include <vs_gen_lib.h>

//-----------------------------------------------------------------------------
// global library internal state variables
vs_gen vs;


////-----------------------------------------------------------------------------
//void init_from_file(std::string filename)
//{
//
//
//}

//-----------------------------------------------------------------------------
void init(/*params*/
    uint32_t sig_tbl_num,
    double* sigma_table_t,
    uint32_t blur_tbl_num,
    uint8_t* dm_values_t,
    uint8_t* br1_table_t,
    uint8_t* br2_table_t,
    uint32_t bg_tbl_num,
    void* bg_tbl_t,
    uint32_t fg_tbl_num,
    void* fg_tbl_t,
    double fg_prob_,
    double bg_prob_,
    uint16_t fg_dm_value_,
    uint16_t bg_dm_value_,
    int32_t max_dm_vals_per_image_
)
{

    vs.init(sig_tbl_num, sigma_table_t, blur_tbl_num, dm_values_t, br1_table_t, br2_table_t, bg_tbl_num, bg_tbl_t, fg_tbl_num,
        fg_tbl_t, fg_prob_, bg_prob_, fg_dm_value_, bg_dm_value_, max_dm_vals_per_image_
    );

}   // end of init


//-----------------------------------------------------------------------------
void generate_scene(unsigned int img_w, 
    unsigned int img_h, 
    unsigned char* img_f1_t, 
    unsigned char* img_f2_t, 
    unsigned char* dm_t
)
{
    uint32_t idx, N;
    int32_t min_N, max_N;
    uint8_t dataset_type = 0;

    std::vector<uint8_t> tmp_br1_table, tmp_br2_table;
    std::vector<uint16_t> dm_vals;
    std::vector<uint16_t> dm_indexes;

    cv::Mat random_img, output_img, mask;
    cv::Mat f1_layer, f2_layer;
    double scale2;

    // assigned the input pointers to cv:Mat containers
    cv::Mat img_f1 = cv::Mat(img_h, img_w, CV_8UC3, img_f1_t);
    cv::Mat img_f2 = cv::Mat(img_h, img_w, CV_8UC3, img_f2_t);
    cv::Mat dm = cv::Mat(img_h, img_w, CV_8UC1, dm_t);
    cv::Size img_size(img_h, img_w);

    // generate random dm_values that include the foreground and background values
    int32_t tmp_dm_num = vs.max_dm_vals_per_image;

    // get the probablility that the background depthmap value will be used
    double bg_x = vs.rng.uniform(0.0, 1.0);

    // get the probability that the foreground depthmap value will be used
    double fg_x = vs.rng.uniform(0.0, 1.0);

    if (bg_x < vs.bg_prob)
        tmp_dm_num--;

    if (fg_x < vs.fg_prob)
        tmp_dm_num--;

    //generate_depthmap_set(min_dm_value, max_dm_value, tmp_dm_num, depthmap_values, dm_values, rng);
    generate_depthmap_index_set(vs.fg_dm_value, vs.bg_dm_value, tmp_dm_num, vs.dm_values, dm_indexes, vs.rng);

    // check the background probability and fill in the tables
    if (bg_x < vs.bg_prob)
    {
        uint16_t dm = vs.rng.uniform(0, vs.bg_br_table.size());
        tmp_br1_table.push_back(vs.bg_br_table[dm].first);
        tmp_br2_table.push_back(vs.bg_br_table[dm].second);

        dm_vals.push_back(vs.bg_dm_value);
    }

    // fill in the tables for the region of interest depthmap values
    for (idx = 0; idx < dm_indexes.size(); ++idx)
    {
        tmp_br1_table.push_back(vs.br1_table[dm_indexes[idx]]);
        tmp_br2_table.push_back(vs.br2_table[dm_indexes[idx]]);
        dm_vals.push_back(vs.dm_values[dm_indexes[idx]]);
    }

    // check the foreground probability and fill in the tables
    if (fg_x < vs.fg_prob)
    {
        uint16_t dm = vs.rng.uniform(0, vs.fg_br_table.size());
        tmp_br1_table.push_back(vs.fg_br_table[dm].first);
        tmp_br2_table.push_back(vs.fg_br_table[dm].second);

        dm_vals.push_back(vs.fg_dm_value);
    }

    //N = (uint32_t)(img_h * img_w * 0.001);
    switch (dataset_type)
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
    dm = cv::Mat(img_h, img_w, CV_8UC1, cv::Scalar::all(dm_vals[0]));

    // blur imgs using dm_values and random masks
    for (idx = 1; idx < dm_vals.size(); ++idx)
    {
        f1_layer = img_f1.clone();
        f2_layer = img_f2.clone();

        // This part help shape the final distribution of depthmap values
        //min_N = (int32_t)ceil((num_objects) / (1 + exp(-0.35 * dm_values[idx] + (0.035 * num_objects))) + 2);
        //min_N = (int32_t)(num_objects / (double)(1.0 + exp(-0.1 * (depthmap_values[dm_indexes[idx]] - (max_dm_value-min_dm_value)/2.0))) + 5);
        //min_N = (int32_t)ceil(((max_dm_value) / (double)(1.0 + exp(-0.35 * depthmap_values[dm_values[idx]] + (0.175 * max_dm_value)))) + 3);
        min_N = (int32_t)ceil(((vs.bg_dm_value) / (double)(1.0 + exp(-0.365 * dm_vals[idx] + (0.175 * vs.bg_dm_value)))) + 3);
        max_N = (int32_t)ceil(2.0 * min_N);  // 2.0

        N = vs.rng.uniform(min_N, max_N + 1);

        // define the scale factor
        scale2 = 60.0 / (double)img_size.width;

        // 
        switch (dataset_type)
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
        overlay_depthmap(dm, mask, dm_vals[idx]);

        // blur f1
        blur_layer(f1_layer, img_f1, mask, vs.blur_kernels[tmp_br1_table[idx]], vs.rng);

        // blur f2
        blur_layer(f2_layer, img_f2, mask, vs.blur_kernels[tmp_br2_table[idx]], vs.rng);
    }

    img_f1_t = img_f1.data;
    img_f2_t = img_f2.data;
    dm_t = dm.data;

}   // end of generate_scene

