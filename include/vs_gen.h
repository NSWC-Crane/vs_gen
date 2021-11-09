#ifndef _VS_GEN_HEADER_H_
#define _VS_GEN_HEADER_H_

// C/C++ includes
#include <cmath>
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
class vs_gen
{
public:

    cv::RNG rng(time(NULL));

    std::vector<uint16_t> dm_values;
    std::vector<uint16_t> dm_indexes;
    int32_t min_N, max_N;

    double scale = 0.1;
    uint32_t BN = 1000;

    std::pair<uint8_t, double> bg_dm;
    std::pair<uint8_t, double> fg_dm;
    std::vector<std::pair<uint8_t, uint8_t>> bg_br_table;
    std::vector<std::pair<uint8_t, uint8_t>> fg_br_table;
    double prob_bg = 0.31;    // set the probablility of selecting the background depthmap value
    double prob_fg = 0.35;    // set the probability of selecting the foreground depthmap value

    //std::vector<cv::Mat> blur_kernels;
    //std::vector<cv::Mat> fft_blur_kernels;

    std::vector<uint8_t> depthmap_values;
    std::vector<double> sigma_table;
    std::vector<uint8_t> br1_table;
    std::vector<uint8_t> br2_table;
    uint8_t dataset_type = 0;

    vs_gen() = default;


    ~vs_gen() 
    {

    }

private:




};




#endif  // _VS_GEN_HEADER_H_
