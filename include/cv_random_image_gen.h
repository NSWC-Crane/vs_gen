#ifndef _CV_RANDOM_IMAGE_GEN_H_
#define _CV_RANDOM_IMAGE_GEN_H_

#include <cstring>

// OpenCV includes
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>


inline void generate_random_shape(cv::Mat& img,
    cv::RNG& rng,
    long nr,
    long nc,
    cv::Scalar color, 
    double scale
)
{
    long x, y;
    long r1, r2, h, w, s;
    double a, angle;
    int32_t radius, max_radius;

    int32_t min_dim = std::min(nr, nc);

    // generate the random point
    x = rng.uniform(0, nc);
    y = rng.uniform(0, nr);

    cv::RotatedRect rect;
    cv::Point2f vertices2f[4];
    cv::Point vertices[4];
    std::vector<cv::Point> pts;
    std::vector<std::vector<cv::Point> > vpts(1);

    // get the shape type
    switch (rng.uniform(0, 3))
    {
    // filled ellipse
    case 0:

        // pick a random radi for the ellipse
        r1 = (long)std::floor(0.5 * scale * rng.uniform(min_dim >> 2, min_dim));
        r2 = (long)std::floor(0.5 * scale * rng.uniform(min_dim >> 2, min_dim));
        a = rng.uniform(0.0, 360.0);

        cv::ellipse(img, cv::Point(x, y), cv::Size(r1, r2), a, 0.0, 360.0, color, -1, cv::LineTypes::LINE_8, 0);
        break;

    // filled rectangle
    case 1:

        h = (long)std::floor(scale * rng.uniform(min_dim >> 2, min_dim));
        w = (long)std::floor(scale * rng.uniform(min_dim >> 2, min_dim));
        a = rng.uniform(0.0, 360.0);

        // Create the rotated rectangle
        rect = cv::RotatedRect(cv::Point(x, y), cv::Size(w, h), (float)a);

        // We take the edges that OpenCV calculated for us
        rect.points(vertices2f);

        // Convert them so we can use them in a fillConvexPoly
        for (int jdx = 0; jdx < 4; ++jdx)
        {
            vertices[jdx] = vertices2f[jdx];
        }

        // Now we can fill the rotated rectangle with our specified color
        cv::fillConvexPoly(img, vertices, 4, color);
        break;

    // 3 to 8 sided filled polygon
    case 2:

        h = (long)std::floor(scale * rng.uniform(min_dim >> 2, min_dim));
        w = (long)std::floor(scale * rng.uniform(min_dim >> 2, min_dim));

        s = rng.uniform(3, 9);
        a = 360.0 / (double)s;

        pts.clear();
        for (long jdx = 0; jdx < s; ++jdx)
        {
            angle = rng.uniform(jdx * a, (jdx + 1) * a);

            if (w / std::abs(std::cos((CV_PI / 180.0) * angle)) <= h / std::abs(std::sin((CV_PI / 180.0) * angle)))
            {
                max_radius = (int32_t)std::abs(w / (double)std::cos((CV_PI / 180.0) * angle));
            }
            else
            {
                max_radius = (int32_t)std::abs(h / (double)std::sin((CV_PI / 180.0) * angle));
            }

            radius = rng.uniform(max_radius >> 2, max_radius);
            pts.push_back(cv::Point((int32_t)(radius * std::cos((CV_PI / 180.0) * angle)), (int32_t)(radius * std::sin((CV_PI / 180.0) * angle))));
        }

        vpts[0] = pts;
        cv::fillPoly(img, vpts, color, cv::LineTypes::LINE_8, 0, cv::Point(x, y));

        break;

    }   // end switch

}   // end of generate_random_shape


// ----------------------------------------------------------------------------------------
void generate_random_image(
    cv::Mat& img,
    cv::RNG& rng,
    long nr, 
    long nc, 
    unsigned int N, 
    double scale   
)
{
    unsigned int idx;

    //long x, y;
    //long r1, r2, h, w, s;
    //double a, angle;    
    //int32_t radius, max_radius;

    //long min_dim = std::min(nr, nc);

    // get the random background color
    cv::Scalar bg_color = cv::Scalar(rng.uniform(0, 256), rng.uniform(0, 256), rng.uniform(0, 256));

    // create the image with the random background color
    img = cv::Mat(nr, nc, CV_8UC3, bg_color);

    // create N shapes
    for (idx = 0; idx < N; ++idx)
    {
        // get the random color for the shape
        cv::Scalar C = cv::Scalar(rng.uniform(0, 256), rng.uniform(0, 256), rng.uniform(0, 256));

        // make sure the color picked is not the background color
        while (C == bg_color)
        {
            C = cv::Scalar(rng.uniform(0, 256), rng.uniform(0, 256), rng.uniform(0, 256));
        }

        generate_random_shape(img, rng, nr, nc, C, scale);

        //// generate the random point
        //x = rng.uniform(0, nc);
        //y = rng.uniform(0, nr);


        //cv::RotatedRect rect;
        //cv::Point2f vertices2f[4];
        //cv::Point vertices[4];
        //std::vector<cv::Point> pts;
        //std::vector<std::vector<cv::Point> > vpts(1);
        ////vpts.push_back(pts);


        ////long x_1; //= -window_width / 2;
        ////long x_2; //= window_width * 3 / 2;
        ////long y_1; //= -window_width / 2;
        ////long y_2; //= window_width * 3 / 2;

        //// get the shape type
        //switch (rng.uniform(0, 3))
        //{
        //// filled ellipse
        //case 0:
        //    
        //    // pick a random radi for the ellipse
        //    r1 = (long)std::floor(0.5 * scale * rng.uniform(min_dim >> 2, min_dim));
        //    r2 = (long)std::floor(0.5 * scale * rng.uniform(min_dim >> 2, min_dim));
        //    a = rng.uniform(0.0, 360.0);

        //    cv::ellipse(img, cv::Point(x, y), cv::Size(r1, r2), a, 0.0, 360.0, C, -1, cv::LineTypes::LINE_8, 0);
        //    break;

        //// filled rectangle
        //case 1:

        //    h = (long)std::floor(scale * rng.uniform(min_dim >> 2, min_dim));
        //    w = (long)std::floor(scale * rng.uniform(min_dim >> 2, min_dim));
        //    a = rng.uniform(0.0, 360.0);

        //    // Create the rotated rectangle
        //    rect = cv::RotatedRect(cv::Point(x,y), cv::Size(w,h), (float)a);

        //    // We take the edges that OpenCV calculated for us
        //    rect.points(vertices2f);

        //    // Convert them so we can use them in a fillConvexPoly
        //    for (int jdx = 0; jdx < 4; ++jdx) 
        //    {
        //        vertices[jdx] = vertices2f[jdx];
        //    }

        //    // Now we can fill the rotated rectangle with our specified color
        //    cv::fillConvexPoly(img, vertices, 4, C);
        //    break;

        //// 3 to 8 sided filled polygon
        //case 2:

        //    h = (long)std::floor(scale * rng.uniform(min_dim >> 2, min_dim));
        //    w = (long)std::floor(scale * rng.uniform(min_dim >> 2, min_dim));

        //    s = rng.uniform(3, 9);
        //    a = 360.0 / (double)s;

        //    cv::Point center(rng.uniform(0, nr), rng.uniform(0, nc));

        //    pts.clear();
        //    for (jdx = 0; jdx < s; ++jdx)
        //    {
        //        angle = rng.uniform(jdx * a, (jdx + 1) * a);

        //        if (w / std::abs(std::cos((CV_PI / 180.0) * angle)) <= h / std::abs(std::sin((CV_PI / 180.0) * angle)))
        //        {
        //            max_radius = (int32_t)std::abs(w / (double)std::cos((CV_PI / 180.0) * angle));
        //        }
        //        else
        //        {
        //            max_radius = (int32_t)std::abs(h / (double)std::sin((CV_PI / 180.0) * angle));
        //        }

        //        radius = rng.uniform(max_radius >> 2, max_radius);
        //        pts.push_back(cv::Point((int32_t)(radius * std::cos((CV_PI / 180.0) * angle)), (int32_t)(radius * std::sin((CV_PI / 180.0) * angle))));
        //    }

        //    vpts[0] = pts;
        //    cv::fillPoly(img, vpts, C, cv::LineTypes::LINE_8, 0);

        //    break;

        //}   // end switch

    }   // end for loop

}   // end of generate_random_image


void generate_random_image(unsigned char*& img,
    long long seed,
    long nr,
    long nc,
    unsigned int N,
    double scale
)
{
    cv::RNG rng(seed);

    cv::Mat cv_img;

    generate_random_image(cv_img, rng, nr, nc, N, scale);

    //    memcpy((void*)data_params, &network_output(0, 0), network_output.size() * sizeof(float));
    img = new unsigned char[cv_img.total()*3];
    std::memcpy((void*)img, cv_img.ptr<unsigned char>(0), cv_img.total()*3);

}   // end of generate_random_image

#endif  // _CV_RANDOM_IMAGE_GEN_H_
