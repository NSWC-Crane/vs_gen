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

#include <turbulence_param.h>
#include <turbulence_sim.h>

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

    std::vector<turbulence_param> roi_tp;
    std::vector<turbulence_param> fg_tp, bg_tp;

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
    std::vector<double> ranges;
    std::vector<double> fg_ranges, bg_ranges;

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
        std::string range_str;

        std::vector<uint8_t> bg_table_br1, bg_table_br2;
        std::vector<uint8_t> fg_table_br1, fg_table_br2;

        try {

            std::ifstream tmp_stream(param_filename);
            std::stringstream buffer;
            buffer << tmp_stream.rdbuf();
            std::string contents = buffer.str();
            tmp_stream.close();

            ryml::Tree config = ryml::parse_in_arena(ryml::to_csubstr(contents));

            // background: depthmap value, probablility, blur radius values
            ryml::NodeRef background = config["background"];
            background["value"] >> bg_dm_value;
            background["probability"] >> bg_prob;
            background["blur_radius1"] >> bg_table_br1;
            background["blur_radius2"] >> bg_table_br2;
            background["ranges"] >> range_str;
            parse_input_range(range_str, bg_ranges);

            vector_to_pair(bg_table_br1, bg_table_br2, bg_br_table);

            // foreground: depthmap value, probablility, blur radius values
            ryml::NodeRef foreground = config["foreground"];
            foreground["value"] >> fg_dm_value;
            foreground["probability"] >> fg_prob;
            foreground["blur_radius1"] >> fg_table_br1;
            foreground["blur_radius2"] >> fg_table_br2;
            foreground["ranges"] >> range_str;
            parse_input_range(range_str, fg_ranges);

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
            roi["ranges"] >> range_str;
            parse_input_range(range_str, ranges);

            // maximum number of depthmap values within a single image
            config["max_dm_values"] >> max_dm_vals_per_image;

            // scaling parameters: pattern scale, shape scale
            ryml::NodeRef scaling_params = config["scaling_params"];
            scaling_params["pattern_scale"] >> pattern_scale;
            scaling_params["shape_scale"] >> shape_scale;

            // final blur sigma index
            config["final_blur_index"] >> final_blur_index;

            ryml::NodeRef turb_params = config["turbulence_parameters"];
            turb_params["aperature"] >> D;

            //ryml::NodeRef range_params = turb_params["ranges"];
            //double min_range, max_range, range_step;
            //range_params["min"] >> min_range;
            //range_params["max"] >> max_range;
            //range_params["step"] >> range_step;

            //generate_range(min_range, max_range, range_step, ranges);

            ryml::NodeRef cn2_params = turb_params["Cn2"];
            cn2_params["min"] >> min_cn2;
            cn2_params["max"] >> max_cn2;

            int bp = 1;

        }
        catch (std::exception &e)
        {
            std::string error_string = "Error parsing input file: " + param_filename + " - " + std::string(e.what()) + "\n";
            error_string += "File: " + std::string(__FILE__) + ", Line #: " + std::to_string(__LINE__);
            throw std::runtime_error(error_string);
        }

        // std::vector<std::vector<std::string>> params;
        // parse_csv_file(param_filename, params);

        auto t = time(NULL);
        rng = cv::RNG(t);

        //std::cout << "seed: " << t << std::endl;

        // precalculate the blur kernels
        generate_blur_kernels();

    }   // end of read_params

    void init_turbulence_params(unsigned int N_)
    {
        uint32_t idx;
        double Cn2 = rng.uniform(min_cn2, max_cn2);
        double obj_size;// = N_ * turbulence_param::get_pixel_size(zoom, L);

        for (idx = 0; idx < fg_ranges.size(); ++idx)
        {
            obj_size = N_* turbulence_param::get_pixel_size(zoom, fg_ranges[idx]);
            fg_tp.push_back(turbulence_param(N_, D, fg_ranges[idx], Cn2, green_wvl, obj_size));
        }
        
        for (idx = 0; idx < bg_ranges.size(); ++idx)
        {
            obj_size = N_ * turbulence_param::get_pixel_size(zoom, bg_ranges[idx]);
            bg_tp.push_back(turbulence_param(N_, D, bg_ranges[idx], Cn2, green_wvl, obj_size));
        }

        for (idx = 0; idx < ranges.size(); ++idx)
        {
            obj_size = N_ * turbulence_param::get_pixel_size(zoom, ranges[idx]);
            roi_tp.push_back(turbulence_param(N_, D, ranges[idx], Cn2, green_wvl, obj_size));
        }
        
    }   // end of init_turbulence_params

//-----------------------------------------------------------------------------
private:

    // gaussian kernel size
    const uint32_t kernel_size = 49;
    double D;
    double min_cn2, max_cn2;
    const uint32_t zoom = 2000;
    const double red_wvl = 525e-9;
    const double green_wvl = 525e-9;
    const double blue_wvl = 525e-9;

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
