#include "Triangle.hpp"
#include "rasterizer.hpp"
#include <cmath>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <chrono>

#define TICK(x) auto bench_##x = std::chrono::steady_clock::now();
#define TOCK(x) std::printf("%s: %lfs\n", #x, std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - bench_##x).count());

// constexpr double MY_PI = 3.1415926;
const float Pi = acos(-1);

Eigen::Matrix4f get_view_matrix(Eigen::Vector3f eye_pos)
{
    Eigen::Matrix4f view = Eigen::Matrix4f::Identity();

    Eigen::Matrix4f translate;
    translate << 1, 0, 0, -eye_pos[0], 0, 1, 0, -eye_pos[1], 0, 0, 1,
        -eye_pos[2], 0, 0, 0, 1;

    view = translate * view;

    return view;
}

Eigen::Matrix4f get_model_matrix(float rotation_angle)
{
    Eigen::Matrix4f model = Eigen::Matrix4f::Identity();
    // TODO: Implement this function
    // Create the model matrix for rotating the triangle around the Z axis.
    // Then return it.

    float rad = rotation_angle / 180.0 * Pi;
    float A = std::cos(rad), B = std::sin(rad);
    model(0, 0) =  A;
    model(0, 1) = -B;
    model(1, 1) =  A;
    model(1, 0) =  B;

    return model;
}


Eigen::Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio,
                                      float zNear, float zFar)
{
    Eigen::Matrix4f projection = Eigen::Matrix4f::Identity();
    // TODO: Implement this function
    // Create the projection matrix for the given parameters.
    // Then return it.

    float n = zNear, f = zFar;
    float fovy = eye_fov / 180.0 * Pi;
    float t = std::abs(n) * std::tan(fovy / 2);
    float b = -t;
    float r = t * aspect_ratio;     // aspect=r/t;
    float l = -r;

    Eigen::Matrix4f ortho_t, ortho_s, persp_to_ortho;
    // persp = ortho * persp_to_ortho
    ortho_s << 2/(r-l), 0, 0, 0,
               0, 2/(t-b), 0, 0,
               0, 0, 2/(n-f), 0,
               0, 0, 0, 1;
    // std::cout << "os:" << ortho_s << std::endl;
    ortho_t << 1, 0, 0, -(r+l)/2,
               0, 1, 0, -(t+b)/2,
               0, 0, 1, -(n+f)/2,
               0, 0, 0, 1;
    // std::cout << "ot:" << ortho_t << std::endl;
    persp_to_ortho << n, 0, 0, 0,
                      0, n, 0, 0,
                      0, 0, n+f, -n*f,
                      0, 0, 1, 0;
    //std::cout << persp_to_ortho << std::endl;
    return ortho_s * ortho_t * persp_to_ortho * projection;
}

Eigen::Matrix4f get_rotation(Eigen::Vector3f axis, float angle) {
    // 提高题：该函数的作用是得到绕任意过原点的轴的旋转变换矩阵。 
    Eigen::Matrix4f rotation = Eigen::Matrix4f::Identity();
    // Rodrigues' Rotation Formula
    float alpha = angle / 180.0 * Pi;

    Matrix3f I = Matrix3f::Identity(), N;
    N << 0, -axis.z(), axis.y(),
         axis.z(), 0, -axis.x(),
         -axis.y(), axis.x(), 0;
    Matrix3f rod = cos(alpha) * I + (1 - cos(alpha)) * axis * axis.transpose() + sin(alpha) * N;
    rotation.block(0, 0, 3, 3) = rod;

    return rotation;
}

int main(int argc, const char** argv)
{
    float angle = 0;
    bool command_line = false;
    std::string filename = "output.png";

    if (argc >= 3) {
        command_line = true;
        angle = std::stof(argv[2]); // -r by default
        if (argc == 4) {
            filename = std::string(argv[3]);
        }
    }

    rst::rasterizer r(700, 700);

    Eigen::Vector3f eye_pos = {0, 0, 5};

    std::vector<Eigen::Vector3f> pos{{2, 0, -2}, {0, 2, -2}, {-2, 0, -2}};

    std::vector<Eigen::Vector3i> ind{{0, 1, 2}};

    auto pos_id = r.load_positions(pos);
    auto ind_id = r.load_indices(ind);

    int key = 0;
    int frame_count = 0;

    if (command_line) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, rst::Primitive::Triangle);
        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);

        cv::imwrite(filename, image);

        return 0;
    }
    // float theta = 0;
    while (key != 27) {
        TICK(frametime)
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);
        // if (key == 'r') {
        //     theta += 10;
        //     std::clog << "rotation applied.\n";
        //     rotation = get_rotation({0, 1, 0}, theta);
        // }
        // r.set_model(get_rotation({0,1,0}, theta) * get_model_matrix(angle));

        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, rst::Primitive::Triangle);

        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::imshow("image", image);
        key = cv::waitKey(10);

        std::cout << "frame count: " << frame_count++ << '\n';

        if (key == 'a') {
            angle += 10;
        }
        else if (key == 'd') {
            angle -= 10;
        }
        TOCK(frametime)
    }

    return 0;
}
