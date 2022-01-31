#include <iostream>
#include <fstream>
#include <omp.h>
#include <complex>
#include <vector>
#include <cmath>

using namespace std;

int compute_pixel(int j, int i, int width, int height, int max_recursion_depth)
{
    // Transform integer pixel coordinates into range of values that suit the mandelbrot set (X=[-2,2], Y=[-1,1])
    double x = ((j - (width/2.0))/width)*4;
    double y = ((i - (height/2.0))/height)*2;
    // Initialize complex numbers c and z for recursion
    complex<double> c = complex<double>(x,y);
    complex<double> z = complex<double>(0,0);
    
    for (int counter = 0; counter < max_recursion_depth; counter++) {
       z = std::pow(z,2) + c;
       // If sequence is not bound (abs(z)>2), then return grayscale value scaled by recursion depth
       if (std::abs(z) > 2) {
           return (ceil(counter/((double) max_recursion_depth - 1)*255));
       }
    }
    
    return 255;
}

int main() {
    int max_recursion_depth = 1000;
    int width = 2000;
    int height = width/2;
    const string image_name = "mandelbrot.pgm";
    vector<string> look_up{256};
    
    for (int i = 0; i < 256; ++i) {
        look_up[i] = to_string(i) + "\n";
    }
    remove(image_name.c_str());
    const double start = omp_get_wtime();
    ofstream image(image_name);
    if (image.is_open()) {
        image << "P2\n" << width << " " << height << " 255\n"; // pgm header
        
#pragma omp parallel
{
        string buffer;
        buffer.reserve(width * 4);
#pragma omp for schedule(dynamic) ordered
        for (int i = 0; i < height; i++) {
            buffer.clear(); 
            for (int j = 0; j < width; j++) {
                i = i;
                j = j;
                buffer += look_up[compute_pixel(j, i, width, height, max_recursion_depth)];
            }
#pragma omp ordered
            image << buffer;
        }
}
        image.close(); // close file output stream
    } else {
        cout << "Could not open the file!";
    }
    cout << omp_get_wtime() - start << " seconds" << endl;
}
