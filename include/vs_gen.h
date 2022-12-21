#ifndef _VS_GEN_HEADER_H_
#define _VS_GEN_HEADER_H_

#include <ryml_all.hpp>

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
//#include <array>
#include <algorithm>
//#include <type_traits>
//#include <list>
#include <set>

// OpenCV includes
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

#include <file_parser.h>

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

    //double scale = 0.1;
    double pattern_scale;
    double shape_scale;

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

    uint32_t final_blur_index = 0;

    vs_gen() = default;

    ~vs_gen()
    {

    }

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
        uint16_t bg_dm_value_,
        int32_t max_dm_vals_per_image_
        ) 
    {

        auto t = time(NULL);
        rng = cv::RNG(t);

        //std::cout << "seed: " << t << std::endl;

        fg_prob = fg_prob_;
        bg_prob = bg_prob_;
        fg_dm_value = fg_dm_value_;
        bg_dm_value = bg_dm_value_;
        max_dm_vals_per_image = max_dm_vals_per_image_;

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

    //-----------------------------------------------------------------------------
    inline void read_params(std::string param_filename)
    {
        //uint32_t idx = 0, jdx = 0;

        std::vector<uint8_t> bg_table_br1, bg_table_br2;
        std::vector<uint8_t> fg_table_br1, fg_table_br2;

        try {

            std::ifstream tmp_stream(param_filename);
            std::stringstream buffer;
            buffer << tmp_stream.rdbuf();
            std::string contents = buffer.str();

            ryml::Tree config = ryml::parse_in_arena(ryml::to_csubstr(contents));

            // background: depthmap value, probablility, blur radius values
            ryml::NodeRef background = config["background"];
            background["value"] >> bg_dm_value;
            background["probability"] >> bg_prob;
            background["blur_radius1"] >> bg_table_br1;
            background["blur_radius2"] >> bg_table_br2;

            vector_to_pair(bg_table_br1, bg_table_br2, bg_br_table);

            // foreground: depthmap value, probablility, blur radius values
            ryml::NodeRef foreground = config["background"];
            foreground["value"] >> fg_dm_value;
            foreground["probability"] >> fg_prob;
            foreground["blur_radius1"] >> fg_table_br1;
            foreground["blur_radius2"] >> fg_table_br2;

            vector_to_pair(fg_table_br1, fg_table_br2, fg_br_table);

            // sigma values, the number of values should be greater than the number of depthmap values
            sigma_table.clear();
            config["sigma"] >> sigma_table;

            // ROI depthmap values, blur radius table 1, blur radius table 2
            dm_values.clear();
            br1_table.clear();
            br2_table.clear();
            ryml::NodeRef roi = config["roi"];
            roi["values"] >> dm_values;
            roi["blur_radius1"] >> br1_table;
            roi["blur_radius2"] >> br2_table;

            // maximum number of depthmap values within a single image
            config["max_dm_values"] >> max_dm_vals_per_image;

            // scaling parameters: pattern scale, shape scale
            ryml::NodeRef scaling_params = config["scaling_params"];
            scaling_params["pattern_scale"] >> pattern_scale;
            scaling_params["shape_scale"] >> shape_scale;

            // final blur sigma index
            config["final_blur_index"] >> final_blur_index;

        }
        catch (std::exception &e)
        {
            throw std::runtime_error("Error parsing input file: " + std::string(e.what()));
        }

        // std::vector<std::vector<std::string>> params;
        // parse_csv_file(param_filename, params);

        auto t = time(NULL);
        rng = cv::RNG(t);

        //std::cout << "seed: " << t << std::endl;

        // precalculate the blur kernels
        generate_blur_kernels();

    }   // end of read_params

//-----------------------------------------------------------------------------
private:

    // gaussian kernel size
    const uint32_t kernel_size = 69;
    //const uint32_t kernel_size = 512;

    //-----------------------------------------------------------------------------
    template<typename T>
    inline void vector_to_pair(std::vector<T>& v1, std::vector<T>& v2, std::vector<std::pair<T, T>>& p1)
    {
        assert(v1.size() == v2.size());

        uint64_t idx;

        p1.clear();

        for (idx = 0; idx < v1.size(); ++idx)
        {
            p1.push_back(std::make_pair(v1[idx], v2[idx]));
        }

    }   // end of vector_to_pair

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
