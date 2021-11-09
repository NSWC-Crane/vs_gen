#ifndef _CV_BLUR_PROCESS_H_
#define _CV_BLUR_PROCESS_H_

#include <cstdint>
#include <set>

// OpenCV includes
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/calib3d.hpp>

#include <cv_random_image_gen.h>



//-----------------------------------------------------------------------------
inline void blur_layer(cv::Mat& input_img,
    cv::Mat &output_img,
    cv::Mat mask,
    cv::Mat &kernel,
    cv::RNG &rng
)
{
    cv::Mat L1_1, L1_2;
    input_img.convertTo(input_img, CV_32FC3);
    output_img.convertTo(output_img, CV_32FC3);
    mask.convertTo(mask, CV_32FC3);

    // blur the src_clone image with the overlay and blur the mask image
    cv::filter2D(input_img, L1_1, -1, kernel, cv::Point(-1, -1), 0.0, cv::BorderTypes::BORDER_REPLICATE);
    cv::filter2D(mask, mask, -1, kernel, cv::Point(-1, -1), 0.0, cv::BorderTypes::BORDER_REPLICATE);

    // multiply the src image times (cv::Scalar(1.0, 1.0, 1.0) - mask)
    cv::multiply(L1_1, mask, L1_1);

    cv::multiply(output_img, cv::Scalar(1.0, 1.0, 1.0) - mask, L1_2);

    // set src equal to L1_1 + L1_2
    cv::add(L1_1, L1_2, output_img);

    output_img.convertTo(output_img, CV_8UC3);

}   // end of blur_layer 


//-----------------------------------------------------------------------------
void generate_depthmap_set(uint16_t min_dm_value, uint16_t max_dm_value, int32_t max_dm_num, 
    const std::vector<uint8_t> depthmap_values, std::vector<uint16_t> &dm_values, cv::RNG rng)
{
    std::set<uint16_t, std::greater<uint16_t>> set_values;

    for (int32_t idx = 0; idx<max_dm_num; idx++)
    {
        uint16_t random_idx = rng.uniform(0, (int32_t)depthmap_values.size());
        set_values.insert(depthmap_values.at(random_idx));
    }

    dm_values = std::vector<uint16_t>(set_values.begin(), set_values.end());
}

//-----------------------------------------------------------------------------
void generate_depthmap_index_set(uint16_t min_dm_value, uint16_t max_dm_value, int32_t max_dm_num,
    const std::vector<uint8_t> depthmap_values, std::vector<uint16_t>& dm_values, cv::RNG rng)
{
    std::set<uint16_t, std::greater<uint16_t>> set_values;

    for (int32_t idx = 0; idx < max_dm_num; idx++)
    {
        uint16_t random_idx = rng.uniform(0, (int32_t)depthmap_values.size());
        set_values.insert(random_idx);
    }

    dm_values = std::vector<uint16_t>(set_values.begin(), set_values.end());
}

//-----------------------------------------------------------------------------
void generate_random_mask(cv::Mat& output_mask,
    cv::Size img_size,
    cv::RNG& rng,
    uint32_t num_shapes,
    double scale = 0.1
    )
{

    uint32_t idx;

    int nr = img_size.width;
    int nc = img_size.height;

    // create the image with a black background color
    output_mask = cv::Mat(nr, nc, CV_8UC3, cv::Scalar::all(0));

    // create N shapes
    for (idx = 0; idx < num_shapes; ++idx)
    {

        // color for all shapes
        cv::Scalar C = cv::Scalar(1, 1, 1);

        // generate the random shape
        generate_random_shape(output_mask, rng, nr, nc, C, scale);

    }   // end for loop

} // end of generate_random_mask


//-----------------------------------------------------------------------------
void generate_random_overlay(cv::Mat random_img, //cv::Size img_size,
    cv::RNG &rng, 
    cv::Mat &output_img, 
    cv::Mat &output_mask, 
    uint32_t num_shapes, 
    double scale = 0.1)
{

    // generate random mask
    generate_random_mask(output_mask, cv::Size(random_img.cols, random_img.rows) , rng, num_shapes, scale*1.3);

    // multiply random_img times output_mask
    cv::multiply(random_img, output_mask, output_img);

} // end of generate_random_overlay


//-----------------------------------------------------------------------------
void overlay_depthmap(cv::Mat &depth_map, cv::Mat mask, uint16_t dm_value)
{
    cv::Mat bg_mask;
    
    cv::cvtColor(mask, mask, cv::COLOR_BGR2GRAY); 

    // set bg_mask = 1 - mask
    cv::subtract(cv::Scalar(1), mask, bg_mask);

    cv::multiply(depth_map, bg_mask, depth_map);

    // multiply mask with dm_value  
    cv::multiply(mask, cv::Scalar(dm_value), mask);

    // overlay depthmap
    cv::add(depth_map, mask, depth_map);

}   // end of overlay_depthmap

//-----------------------------------------------------------------------------
void overlay_image(cv::Mat& input_img, cv::Mat &overlay, cv::Mat mask)
{

    // multiply input image by 1-mask, mask should be 0 or 1, 1's should be where the overlay goes
    cv::multiply(input_img, cv::Scalar(1, 1, 1) - mask, input_img);

    // overlay the image.  Assuming that the overlay has already been zeroed in the right spots
    cv::add(input_img, overlay, input_img);

}   // end of overlay_image

//-----------------------------------------------------------------------------
void new_shapes(cv::Mat &img, uint32_t img_h, uint32_t img_w, cv::RNG rng)
{
    img = cv::Mat(img_h, img_w, CV_8UC1, cv::Scalar::all(0));
    std::vector<cv::Point> pts;
    std::vector<std::vector<cv::Point> > vpts(1);
    double scale = 0.2;
    int num_shapes = 10;
    long h, w, s;
    double a;
    int radius, max_radius;
    int angle;

    for (int N = 0; N < num_shapes; N++)
    {
        h = (long)std::floor(0.5 * scale * rng.uniform(0, std::min(img_h, img_w)));
        w = (long)std::floor(0.5 * scale * rng.uniform(0, std::min(img_h, img_w)));
        
        s = rng.uniform(3, 9);
        a = 360 / (double)s;
        
        cv::Point center(rng.uniform(0, img_h), rng.uniform(0, img_w));

        pts.clear();
        for (int i = 0; i < s; i++)
        {
            angle = (int32_t)rng.uniform((double)(i * a), (double)((i + 1) * a));
            
            if (w/std::abs(std::cos((CV_PI / 180.0) * angle)) <= h/std::abs(std::sin((CV_PI / 180.0) * angle)))
            {
                max_radius = (int32_t)std::abs(w / (double)std::cos((CV_PI / 180.0) * angle));
            }
            else
            {
                max_radius = (int32_t)std::abs(h / (double)std::sin((CV_PI / 180.0) * angle));
            }
            radius = rng.uniform(max_radius/4, max_radius);
            pts.push_back(cv::Point((int32_t)(radius * std::cos((CV_PI / 180.0) * angle)), (int32_t)(radius * std::sin((CV_PI / 180.0) * angle))));
        }

        vpts[0] = pts;
        cv::fillPoly(img, vpts, cv::Scalar(255), cv::LineTypes::LINE_8, 0, center);
        // display bounded box 
        cv::rectangle(img, cv::Point(center.x - w, center.y - h), cv::Point(center.x + w, center.y + h), cv::Scalar(128));
    }
}

#endif // _CV_BLUR_PROCESS_H_
