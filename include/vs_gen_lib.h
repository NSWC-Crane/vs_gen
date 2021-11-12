#ifndef _VS_GEN_LIB_H_
#define _VS_GEN_LIB_H_

#if defined(_WIN32) | defined(__WIN32__) | defined(__WIN32) | defined(_WIN64) | defined(__WIN64)

    #ifdef LIB_EXPORTS
        #define VS_GEN_LIB __declspec(dllexport)
    #else
        #define VS_GEN_LIB 
//__declspec(dllimport)
    #endif

#else
    #define VS_GEN_LIB

#endif


// ----------------------------------------------------------------------------------------
#ifdef __cplusplus
extern "C" {
#endif
    VS_GEN_LIB void init(uint32_t sig_tbl_num,
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
    );
#ifdef __cplusplus
}
#endif

// ----------------------------------------------------------------------------------------
#ifdef __cplusplus
extern "C" {
#endif
    VS_GEN_LIB void generate_scene(unsigned int img_w, 
        unsigned int img_h, 
        unsigned char* img_f1_t, 
        unsigned char* img_f2_t, 
        unsigned char* dm_t
    );
#ifdef __cplusplus
}
#endif


#endif  // _VS_GEN_LIB_H_
