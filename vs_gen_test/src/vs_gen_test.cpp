#define _CRT_SECURE_NO_WARNINGS

#if defined(_WIN32) | defined(__WIN32__) | defined(__WIN32) | defined(_WIN64) | defined(__WIN64)
#include <windows.h>

#else
#include <dlfcn.h>
typedef void* HINSTANCE;

#endif

// C/C++ includes
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <string>
#include <vector>
#include <algorithm>

// OpenCV includes
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

// custom includes
#include <num2string.h>

#define USE_LIB

#if defined(USE_LIB)
//#include <vs_gen_lib.h>

typedef void (*lib_init_vs_gen_from_file)(const char* fn);
typedef void (*lib_set_vs_seed)(int seed);
typedef void (*lib_generate_vs_scene)(unsigned int img_w, unsigned int img_h, unsigned char* img_f1_t, unsigned char* img_f2_t, unsigned char* dm_t);

#else
#include <vs_gen.h>

#endif

// ----------------------------------------------------------------------------------------
inline std::ostream& operator<<(std::ostream& out, std::vector<uint8_t>& item)
{
    for (uint64_t idx = 0; idx < item.size() - 1; ++idx)
    {
        out << static_cast<uint32_t>(item[idx]) << ",";
    }
    out << static_cast<uint32_t>(item[item.size() - 1]);
    return out;
}

// ----------------------------------------------------------------------------------------
template <typename T>
inline std::ostream& operator<<(std::ostream& out, std::vector<T>& item)
{
    for (uint64_t idx = 0; idx < item.size() - 1; ++idx)
    {
        out << item[idx] << ",";
    }
    out << item[item.size() - 1];
    return out;
}

//-----------------------------------------------------------------------------
int main(int argc, char** argv)
{
    int bp = 0;

    uint32_t idx = 0, jdx = 0;
    uint32_t img_h = 256;
    uint32_t img_w = 256;
    cv::Size img_size(img_h, img_w);

    cv::RNG rng(time(NULL));

    // timing variables
    typedef std::chrono::duration<double> d_sec;
    auto start_time = std::chrono::system_clock::now();
    auto stop_time = std::chrono::system_clock::now();
    auto elapsed_time = std::chrono::duration_cast<d_sec>(stop_time - start_time);

    cv::Mat img_f1, img_f2;
    cv::Mat dm_img;
    cv::Mat montage;

    uint64_t vs_seed;

    std::string window_name = "image";

    std::string lib_filename;

    std::string image_tb_num, image_num;
    std::string fp1_image, fp2_image, dm_image;
    uint32_t num_tb_images = 20;
    uint32_t num_images = 450;

    std::ofstream data_log_stream;
    std::string log_filename = "../results/tb23a_input.txt";


    // setup the windows to display the results
    //cv::namedWindow(window_name, cv::WINDOW_NORMAL);
    //cv::resizeWindow(window_name, 2*img_w, img_h);

    std::vector<int32_t> compression_params;
    compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
    compression_params.push_back(2);

    // do work here
    try
    {    

        // std::cout << std::string(argv[0]) << std::endl;

#if defined(USE_LIB)
    // load in the library
    #if defined(_WIN32) | defined(__WIN32__) | defined(__WIN32) | defined(_WIN64) | defined(__WIN64)
#if defined(_DEBUG)
        lib_filename = "../../vs_gen_lib/build/Debug/vs_gen.dll";
#else
        lib_filename = "../../../vs_gen_lib/build/Release/vs_gen.dll";
#endif

        //HINSTANCE vs_gen_lib = LoadLibraryEx(lib_filename.c_str(), 0, DONT_RESOLVE_DLL_REFERENCES);
        HINSTANCE vs_gen_lib = LoadLibrary(lib_filename.c_str());

        if (vs_gen_lib == NULL)
        {
            int load_error = GetLastError();
            throw std::runtime_error("error loading library - code: " + std::to_string(load_error) + " File: " + std::string(__FILE__) + ", Line #: " + std::to_string(__LINE__));
        }

        lib_init_vs_gen_from_file init_vs_gen_from_file = (lib_init_vs_gen_from_file)GetProcAddress(vs_gen_lib, "init_vs_gen_from_file");
        lib_set_vs_seed set_vs_seed = (lib_set_vs_seed)GetProcAddress(vs_gen_lib, "set_vs_seed");
        lib_generate_vs_scene generate_vs_scene = (lib_generate_vs_scene)GetProcAddress(vs_gen_lib, "generate_vs_scene");

    #else
        lib_filename = "../../vs_gen_lib/build/libvs_gen.so";
        void* vs_gen_lib = dlopen(lib_filename.c_str(), RTLD_NOW);

        if (vs_gen_lib == NULL)
        {
            throw std::runtime_error("error loading library");
        }

        //lib_init_turbulence_params init_turbulence_params = (lib_init_turbulence_params)dlsym(turb_lib, "init_turbulence_params");
        //lib_apply_turbulence apply_turbulence = (lib_apply_turbulence)dlsym(turb_lib, "apply_turbulence");


        lib_init_vs_gen_from_file init_vs_gen_from_file = (lib_init_vs_gen_from_file)dlsym(vs_gen_lib, "init_vs_gen_from_file");
        lib_set_vs_seed set_vs_seed = (lib_set_vs_seed)dlsym(vs_gen_lib, "set_vs_seed");
        lib_generate_vs_scene generate_vs_scene = (lib_generate_vs_scene)dlsym(vs_gen_lib, "generate_vs_scene");

    #endif

#endif
        std::string input_filename = "../../blur_params_v23a.yml";

#if defined(USE_LIB)
        init_vs_gen_from_file(input_filename.c_str());
#else

#endif
        data_log_stream.open((log_filename), std::ios::out);

        // Add the date and time to the start of the log file
        data_log_stream << std::endl << "#------------------------------------------------------------------------------" << std::endl;

        //-----------------------------------------------------------------------------
        cv::Mat img_tilt;
        img_f1 = cv::Mat::zeros(img_h, img_w, CV_8UC3);
        img_f2 = cv::Mat::zeros(img_h, img_w, CV_8UC3);
        dm_img = cv::Mat::zeros(img_h, img_w, CV_8UC1);

        char key = 0;
        //cv::resizeWindow(window_name, 4*N, 2*N);


        for (jdx = 0; jdx <= num_images; ++jdx)
        {
            vs_seed = time(NULL);

            image_num = num2str(jdx, "%04d");

            std::cout << "image #" << image_num << "..." << std::endl;

            //while(key != 'q')
            for (idx = 0; idx < num_tb_images; ++idx)
            {

                set_vs_seed(vs_seed);
                start_time = std::chrono::system_clock::now();

#if defined(USE_LIB)

                // apply_turbulence(N, N, img.ptr<double>(0), img_blur.ptr<double>(0));
                generate_vs_scene(img_w, img_h, img_f1.ptr<uint8_t>(0), img_f2.ptr<uint8_t>(0), dm_img.ptr<uint8_t>(0));

#else
                //generate_tilt_image(img, Pv[0], rng, img_tilt);
                //generate_blur_image(img_tilt, Pv[0], rng, img_blur);

                //generate_tilt_image(img, Pv[22], rng, img_tilt);
                //generate_blur_image(img_tilt, Pv[22], rng, img_blur2);
#endif

            //img_blur.convertTo(img_blur, CV_8UC1);

                stop_time = std::chrono::system_clock::now();
                elapsed_time = std::chrono::duration_cast<d_sec>(stop_time - start_time);

                image_tb_num = num2str(idx, "%02d");
                std::cout << "tb #" << image_tb_num + " time (s): " << elapsed_time.count() << std::endl;

                //cv::hconcat(img_f1, img_f2, montage);
                //cv::hconcat(montage, dm_img, montage);
                //cv::imshow(window_name, montage);
                //cv::imshow("Depth Map", dm_img*10);
                //key = cv::waitKey(0);

                fp1_image = "images/test_image_fp1_i" + image_num + "_tb" + image_tb_num + ".png";
                fp2_image = "images/test_image_fp2_i" + image_num + "_tb" + image_tb_num + ".png";
                dm_image = "depth_maps/test_image_dm_i" + image_num + ".png";

                cv::imwrite("../results/" + fp1_image, img_f1, compression_params);
                cv::imwrite("../results/" + fp2_image, img_f2, compression_params);

                data_log_stream << fp1_image + ", " + fp2_image + ", " + dm_image << std::endl;
            }   // end idx

            cv::imwrite("../results/" + dm_image, dm_img, compression_params);
        
        }   // end jdx


    }
    catch(std::exception& e)
    {
        std::cout << "Error: " << e.what() << std::endl;
        std::cout << "Filename: " << __FILE__ << ", Line #: " << __LINE__ << std::endl;
    }

    std::cout << std::endl << "End of Program.  Press Enter to close..." << std::endl;
	std::cin.ignore();
    cv::destroyAllWindows();

}   // end of main

