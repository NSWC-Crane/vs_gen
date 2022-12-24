#ifndef _VS_GEN_LIB_H_
#define _VS_GEN_LIB_H_

#if defined(_WIN32) | defined(__WIN32__) | defined(__WIN32) | defined(_WIN64) | defined(__WIN64)

    #ifdef LIB_EXPORTS
        #define VS_GEN_LIB __declspec(dllexport)
    #else
        #define VS_GEN_LIB __declspec(dllimport)
    #endif

#else
    #define VS_GEN_LIB

#endif

//-----------------------------------------------------------------------------
#ifdef __cplusplus
extern "C" {
#endif
    //-----------------------------------------------------------------------------
    VS_GEN_LIB void init_vs_gen_from_file(const char* fn);

    //-----------------------------------------------------------------------------
    VS_GEN_LIB void init_vs_generator(unsigned int sig_tbl_num,
        double* sigma_table_t,
        unsigned int blur_tbl_num,
        unsigned char* dm_values_t,
        unsigned char* br1_table_t,
        unsigned char* br2_table_t,
        unsigned int bg_tbl_num,
        void* bg_tbl_t,
        unsigned int fg_tbl_num,
        void* fg_tbl_t,
        double fg_prob_,
        double bg_prob_,
        unsigned short fg_dm_value_,
        unsigned short bg_dm_value_,
        int max_dm_vals_per_image_
    );

    //-----------------------------------------------------------------------------
    VS_GEN_LIB void set_vs_seed(int seed);

    //-----------------------------------------------------------------------------
    VS_GEN_LIB void get_vs_minmax(unsigned short* min_dm_value, unsigned short* max_dm_value);

    //-----------------------------------------------------------------------------
    VS_GEN_LIB void set_vs_shape_scale(double s);

    //-----------------------------------------------------------------------------
    VS_GEN_LIB double get_vs_shape_scale();

    //-----------------------------------------------------------------------------
    VS_GEN_LIB double get_vs_pattern_scale();

    //-----------------------------------------------------------------------------
    VS_GEN_LIB void set_vs_pattern_scale(double s);

    //-----------------------------------------------------------------------------
    VS_GEN_LIB void generate_vs_scene(
        unsigned int img_w,
        unsigned int img_h, 
        unsigned char* img_f1_t, 
        unsigned char* img_f2_t, 
        unsigned char* dm_t
    );

#ifdef __cplusplus
}
#endif


#endif  // _VS_GEN_LIB_H_
