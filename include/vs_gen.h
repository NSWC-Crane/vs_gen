#ifndef _VS_GEN_HEADER_H_
#define _VS_GEN_HEADER_H_

// C/C++ includes
#include <cmath>
#include <ctime>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
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

//-----------------------------------------------------------------------------
bool compare(std::pair<uint8_t, uint8_t> p1, std::pair<uint8_t, uint8_t> p2)
{
    return std::max(p1.first, p1.second) < std::max(p2.first, p2.second);
}

//-----------------------------------------------------------------------------
class vs_gen
{
public:

    cv::RNG rng;

    std::vector<uint16_t> dm_values;
    std::vector<uint16_t> dm_indexes;
    int32_t min_N, max_N;

    double scale = 0.1;
    uint32_t BN = 1000;

    std::pair<uint16_t, double> bg_dm;
    std::pair<uint16_t, double> fg_dm;
    std::vector<std::pair<uint8_t, uint8_t>> bg_br_table;
    std::vector<std::pair<uint8_t, uint8_t>> fg_br_table;
    double prob_bg = 0.31;      // set the probablility of selecting the background depthmap value
    double prob_fg = 0.35;      // set the probability of selecting the foreground depthmap value
    uint16_t min_dm_value;      // = fg_dm.first;        // depthmap_values.front();
    uint16_t max_dm_value;      // = bg_dm.first;        // depthmap_values.back();

    std::vector<cv::Mat> blur_kernels;
    //std::vector<cv::Mat> fft_blur_kernels;

    std::vector<uint8_t> depthmap_values;
    std::vector<double> sigma_table;
    std::vector<uint8_t> br1_table;
    std::vector<uint8_t> br2_table;
    uint8_t dataset_type = 0;

    int32_t max_dm_vals_per_image;


    vs_gen()
    {
        rng = cv::RNG(time(NULL));









        generate_blur_kernels();

    }



    ~vs_gen() 
    {

    }

private:

    // gaussian kernel size
    const uint32_t kernel_size = 69;
    //const uint32_t kernel_size = 512;

    void generate_blur_kernels(void)
    {    
        uint32_t idx;
        cv::Mat kernel;

        // get the max value of blur radius
        uint8_t max_br_value = *std::max_element(br1_table.begin(), br1_table.end());
        max_br_value = std::max(max_br_value, *std::max_element(br2_table.begin(), br2_table.end()));
        auto p1 = *std::max_element(bg_br_table.begin(), bg_br_table.end(), compare);
        auto p2 = *std::max_element(fg_br_table.begin(), fg_br_table.end(), compare);
        max_br_value = std::max(max_br_value, std::max(p1.first, p1.second));
        max_br_value = std::max(max_br_value, std::max(p2.first, p2.second));

        // gerenate all of the kernels once and then use them throughout
        for (idx = 0; idx <= max_br_value; ++idx)
        {
            create_gaussian_kernel(kernel_size, sigma_table[idx], kernel);
            blur_kernels.push_back(kernel.clone());

            //cv::dft(kernel, tmp);// , cv::DFT_COMPLEX_OUTPUT);
            //fft_blur_kernels.push_back(tmp.clone());
        }
    }   // end of generate_blur_kernels

};




#endif  // _VS_GEN_HEADER_H_
