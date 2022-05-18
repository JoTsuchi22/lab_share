#include <matplotlib-cpp/matplotlibcpp.h>

namespace plt = matplotlibcpp;

int main(){
    int n = 1000;
    std::vector<double> x(n), y(n);

    for(int i=0; i<n; ++i) {
        x[i] = 2 * M_PI * i / n;
        y[i] = sin(x[i]);
    }

    plt::plot(x, y);
    plt::show();
}