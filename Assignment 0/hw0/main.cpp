#include <cmath>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <iostream>


// #define PA0_DEMO

void Pa0Func() {
    /* 
    * PA 0
    */
    // DONE: Define point P
    // DONE: Define rotation matrix M
    // DONE: M * P
    const float PI = 3.1415926f;
    // 2D homo coord
    Eigen::Vector3f p = {2, 1, 1};
    std::cout << "point p^T at:" << p.transpose() << '\n';



    Eigen::Matrix3f m, n;
    // m: rotation matrix
    float t = 45.0 / 180.0 * PI;
    float a = std::sin(t), b = std::cos(t);

    m << b, -a, 0,
         a,  b, 0,
         0,  0, 1;

    // n: translation matrix
    float tx = 1.0f, ty = 2.0f;

    n << 1, 0, tx,
         0, 1, ty,
         0, 0,  1;

    std::cout << "rotation mat, m: \n" << m << '\n';

    std::cout << "m*p=\n" << m * p << '\n';
    auto m_inv = m.inverse();
    std::cout << "m_inv*m*p=\n" << m_inv * m * p << '\n';

    std::cout << "translate only for n * p = \n" << n * p << '\n';

    std::cout << "rotate then translate: n*m*p = \n" << n * m * p << '\n';
}


int main(){

#ifdef PA0_DEMO
    // Basic Example of cpp
    std::cout << "Example of cpp \n";
    float a = 1.0, b = 2.0;
    std::cout << a << std::endl;
    std::cout << a/b << std::endl;
    std::cout << std::sqrt(b) << std::endl;
    std::cout << std::acos(-1) << std::endl;
    std::cout << std::sin(30.0/180.0*acos(-1)) << std::endl;

    // Example of vector
    std::cout << "Example of vector \n";
    // vector definition
    Eigen::Vector3f v(1.0f,2.0f,3.0f);
    Eigen::Vector3f w(1.0f,0.0f,0.0f);
    // vector output
    std::cout << "Example of output \n";
    std::cout << v << std::endl;
    // vector add
    std::cout << "Example of add \n";
    std::cout << v + w << std::endl;
    // vector scalar multiply
    std::cout << "Example of scalar multiply \n";
    std::cout << v * 3.0f << std::endl;
    std::cout << 2.0f * v << std::endl;

    // Example of matrix
    std::cout << "Example of matrix \n";
    // matrix definition
    Eigen::Matrix3f i,j;
    i << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
    j << 2.0, 3.0, 1.0, 4.0, 6.0, 5.0, 9.0, 7.0, 8.0;
    // matrix output
    std::cout << "Example of output \n";
    std::cout << i << std::endl;
    // matrix add i + j
    // matrix scalar multiply i * 2.0
    // matrix multiply i * j
    // matrix multiply vector i * v
#endif
    Pa0Func();

    return 0;
}