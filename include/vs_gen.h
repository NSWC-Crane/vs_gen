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

    //std::vector<uint16_t> dm_values;
    //std::vector<uint16_t> dm_indexes;
    //int32_t min_N, max_N;

    double scale = 0.1;
    uint32_t BN = 1000;

    //std::pair<uint16_t, double> bg_dm;
    //std::pair<uint16_t, double> fg_dm;

    double fg_prob = 0.35;      // the probability of selecting the foreground depthmap value (0.35)
    double bg_prob = 0.31;      // the probablility of selecting the background depthmap value (0.31)
    uint16_t fg_dm_value;      // = fg_dm.first;        // foreground depthmap value; it is assumed that this value is the smallest depthmap value in the set
    uint16_t bg_dm_value;      // = bg_dm.first;        // background depthmap value; it is assumed that this value is the largest depthmap value in the set

    std::vector<cv::Mat> blur_kernels;
    //std::vector<cv::Mat> fft_blur_kernels;

    std::vector<double> sigma_table;
    std::vector<uint8_t> dm_values;
    std::vector<uint8_t> br1_table;
    std::vector<uint8_t> br2_table;
    std::vector<std::pair<uint8_t, uint8_t>> bg_br_table;
    std::vector<std::pair<uint8_t, uint8_t>> fg_br_table;

    int32_t max_dm_vals_per_image;


    vs_gen() = default;

    //vs_gen(
    //    std::vector<double> st_,
    //    std::vector<uint8_t> dm_values_,
    //    std::vector<uint8_t> br1_table_,
    //    std::vector<uint8_t> br2_table_,
    //    std::vector<std::pair<uint8_t, uint8_t>> bg_br_table_,
    //    std::vector<std::pair<uint8_t, uint8_t>> fg_br_table_,
    //    double fg_prob_,
    //    double bg_prob_,
    //    uint16_t fg_dm_value_,
    //    uint16_t bg_dm_value_
    //) : sigma_table(st_), dm_values(dm_values_), br1_table(br1_table_), br2_table(br2_table_), 
    //    bg_br_table(bg_br_table_), fg_br_table(fg_br_table_), fg_prob(fg_prob_), bg_prob(bg_prob_), 
    //    fg_dm_value(fg_dm_value_), bg_dm_value(bg_dm_value_)
    //{
    //    rng = cv::RNG(time(NULL));

    //    generate_blur_kernels();
    //}

    //-----------------------------------------------------------------------------
    inline void init(uint32_t sig_tbl_num,
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
        uint16_t bg_dm_value_
        ) 
    {

        rng = cv::RNG(time(NULL));

        fg_prob = fg_prob_;
        bg_prob = bg_prob_;
        fg_dm_value = fg_dm_value_;
        bg_dm_value = bg_dm_value_;

        // copy the pointer data into the sigma table vector
        sigma_table.insert(sigma_table.end(), &sigma_table_t[0], &sigma_table_t[sig_tbl_num]);

        // copy the pointer data into the depthmap values
        dm_values.insert(dm_values.end(), &dm_values_t[0], &dm_values_t[blur_tbl_num]);

        // copy the pointer data into the blur rate tables
        br1_table.insert(br1_table.end(), &br1_table_t[0], &br1_table_t[blur_tbl_num]);
        br2_table.insert(br2_table.end(), &br2_table_t[0], &br2_table_t[blur_tbl_num]);

        auto bg = (std::pair<uint8_t, uint8_t>*)bg_tbl_t;
        auto fg = (std::pair<uint8_t, uint8_t>*)fg_tbl_t;

        bg_br_table.insert(bg_br_table.end(), &bg[0], &bg[bg_tbl_num]);
        fg_br_table.insert(fg_br_table.end(), &fg[0], &fg[fg_tbl_num]);

        // precalculate the blur kernels
        generate_blur_kernels();

    }

    //vs_gen(vs_gen v_)
    //{

    //}


    ~vs_gen() 
    {

    }

//-----------------------------------------------------------------------------
private:

    // gaussian kernel size
    const uint32_t kernel_size = 69;
    //const uint32_t kernel_size = 512;

    //-----------------------------------------------------------------------------
    // this function precomputes the blur kernels based on the input sigma table and the kerrnel size
    void generate_blur_kernels(void)
    {    
        uint32_t idx;
        cv::Mat kernel;

        // get the max value of blur radius - assume that the sigma table is 0 index based
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
