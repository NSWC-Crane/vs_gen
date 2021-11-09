#ifndef _CV_DFT_CONV_H_
#define _CV_DFT_CONV_H_

#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>


void fftshift(cv::Mat &img)
{
    // get the center
    int cx = img.cols >> 1;
    int cy = img.rows >> 1;

    cv::Mat q0(img, cv::Rect(0, 0, cx, cy));   // Top-Left - Create a ROI per quadrant
    cv::Mat q1(img, cv::Rect(cx, 0, cx, cy));  // Top-Right
    cv::Mat q2(img, cv::Rect(0, cy, cx, cy));  // Bottom-Left
    cv::Mat q3(img, cv::Rect(cx, cy, cx, cy)); // Bottom-Right

    cv::Mat tmp;                            // swap quadrants (Top-Left with Bottom-Right)
    q0.copyTo(tmp);
    q3.copyTo(q0);
    tmp.copyTo(q3);

    q1.copyTo(tmp);                     // swap quadrant (Top-Right with Bottom-Left)
    q2.copyTo(q1);
    tmp.copyTo(q2);
}


void dft_conv_rgb(cv::Mat &img, cv::Mat &dft_kernel, cv::Mat &dst)
{
    uint8_t idx = 0;

    // dft_kernel should already be in DFT form

    // allocate temporary buffers and initialize them with 0's
    cv::Mat tmp;// = cv::Mat::zeros(img.size(), CV_32FC2);

    // split image into the three channels
    std::vector<cv::Mat> rgb(3);
    cv::split(img, rgb);

    for (idx = 0; idx < rgb.size(); ++idx)
    {
        rgb[idx].convertTo(rgb[idx], CV_64FC1);
        // now transform the input image on each color channel
        cv::dft(rgb[idx], tmp); // cv::DFT_COMPLEX_OUTPUT

        // multiply the spectrums; the function handles packed spectrum representations well
        cv::mulSpectrums(tmp, dft_kernel, tmp, 0);

        // transform the product back from the frequency domain.
        cv::dft(tmp, rgb[idx], cv::DFT_INVERSE + cv::DFT_SCALE);

        fftshift(rgb[idx]);
    }

    cv::merge(rgb, dst);

    dst.convertTo(dst, CV_8UC3);

}   // end of dft_conv_rgb


#endif	// _CV_DFT_CONV_H_