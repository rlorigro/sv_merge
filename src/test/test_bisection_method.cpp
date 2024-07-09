#include <functional>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <map>
#include <unordered_map>
#include <utility>
#include <random>

using std::function;
using std::runtime_error;
using std::cerr;
using std::floor;
using std::map;
using std::unordered_map;
using std::pair;
using std::mt19937;
using std::uniform_real_distribution;



double get_normalized_distance(double d, double n, double n_min, double n_max, double d_min, double d_max){
    double n_range = n_max - n_min;
    double d_range = d_max - d_min;

    if (n_range == 0){
        n_range = 1;
    }

    if (d_range == 0){
        d_range = 1;
    }

    // Compute normalized distance from utopia point as "cost"
    double d_norm = (d - d_min)/(d_range);
    double n_norm = (n - n_min)/(n_range);

    return sqrt(d_norm*d_norm + n_norm*n_norm);
}


pair<int64_t,double> solve_with_binary_search(
        function<double(double)> f,
        int64_t n_min,
        int64_t n_max,
        double d_min,
        double d_max
){
    double n_cost;
    int64_t n;

    // With the range of n determined, we can now use binary search to find the optimal n and d values
    unordered_map<int64_t, double> results;

    // Add the extreme values to the results map
    results.emplace(n_min, get_normalized_distance(d_max, n_min, n_min, n_max, d_min, d_max));
    results.emplace(n_max, get_normalized_distance(d_min, n_max, n_min, n_max, d_min, d_max));

    // Set arbitrary limit on maximum iterations
    int64_t max_iter = 20;

    int64_t i = 0;
    auto a_i = n_min;
    auto b_i = n_max;

    double phi_inverse = (sqrt(5) - 1)/2;

    // Minimize d for each given n until it can be proven that the d value is optimal (left and right values are larger)
    while (i < max_iter){
        ///        c = b - (b - a) * invphi
        ///        d = a + (b - a) * invphi
        ///        if f(c) < f(d):
        ///            b = d
        ///        else:  # f(c) > f(d) to find the maximum
        ///            a = c

        int64_t c_i = floor(double(b_i) - (double(b_i - a_i) * phi_inverse));
        int64_t d_i = floor(double(a_i) + (double(b_i - a_i) * phi_inverse));

        double c_cost = -1;
        double d_cost = -1;

        auto c_result = results.find(c_i);
        if (c_result == results.end()){
            double c = f(c_i);
            c_cost = get_normalized_distance(c, c_i, d_min, d_max, n_min, n_max);
            results.emplace(c_i, c_cost);
        }
        else{
            c_cost = c_result->second;
        }

        auto d_result = results.find(d_i);
        if (d_result == results.end()){
            double d = f(d_i);
            d_cost = get_normalized_distance(d, d_i, d_min, d_max, n_min, n_max);
            results.emplace(d_i, d_cost);

        }
        else{
            d_cost = d_result->second;
        }

        if (c_cost < d_cost){
            // Prevent recomputing the same interval
            if (b_i == d_i){
                break;
            }

            b_i = d_i;
        } else {
            // Prevent recomputing the same interval
            if (a_i == c_i){
                break;
            }

            a_i = c_i;
        }

        i++;
    }

    int64_t c = floor(double(b_i) - (double(b_i - a_i) * phi_inverse));
    int64_t d = floor(double(a_i) + (double(b_i - a_i) * phi_inverse));

    double c_cost = results.at(c);
    double d_cost = results.at(d);

    // Compute the final minimum
    if (c_cost < d_cost){
        n = c;
        n_cost = c_cost;
    } else {
        n = d;
        n_cost = d_cost;
    }

    return {n, n_cost};
}


class F{
public:
    double slope;
    double intercept;

    F(double slope, double intercept) : slope(slope), intercept(intercept) {}

    double operator()(double x){
        return slope/(x + intercept);
    }
};


void test_solve(double slope, double intercept, int64_t n_min, int64_t n_max){
    F f(slope, intercept);

    double d_min = f(n_max);
    double d_max = f(n_min);

    double cost_min_empirical = std::numeric_limits<double>::max();
    int64_t n_min_empirical = -1;

    for (auto i = n_min; i <= n_max; i++){
        auto y = f(i);

        // Compute euclidean distance of i and f(i) from utopia point
        double cost = get_normalized_distance(y, i, d_min, d_max, n_min, n_max);

        if (cost < cost_min_empirical){
            cost_min_empirical = cost;
            n_min_empirical = i;
        }
    }

    auto [n,n_cost] = solve_with_binary_search(f, n_min, n_max, d_min, d_max);

    if (n != n_min_empirical){
        // Print all the integer values from 0 to 10
        for (auto i = n_min; i <= n_max; i++){
            auto y = f(i);

            // Compute euclidean distance of i and f(i) from utopia point
            double cost = get_normalized_distance(y, i, d_min, d_max, n_min, n_max);

            cerr << "i: " << i << "\tcost: " << cost << '\n';
        }

        cerr << "n_min_empirical: " << n_min_empirical << '\n';
        cerr << "n_min: " << n << '\n';
        throw runtime_error("Mismatch between empirical and result minimum");
    }
}


int main() {
    // create 4 RNGs for slope, intercept, min, and max
    mt19937 slope_generator(37);
    mt19937 intercept_generator(69);
    mt19937 n_min_generator(37);
    mt19937 n_max_generator(73);

    for (int i = 0; i < 100'000; i++){
        uniform_real_distribution<double> slope_dist(0.0, 2.0);
        uniform_real_distribution<double> intercept_dist(0, 10.0);
        uniform_real_distribution<double> n_min_dist(0, 11);
        uniform_real_distribution<double> n_max_dist(10, 20);


        double slope = slope_dist(slope_generator);
        double intercept = intercept_dist(intercept_generator);
        int64_t n_min = n_min_dist(n_min_generator);
        int64_t n_max = n_max_dist(n_max_generator);

        cerr << "slope:\t" << round(slope*1000)/1000 << "\tintercept:\t" << round(intercept*1000)/1000 << "\tn_min:\t" << n_min << "\tn_max:\t" << n_max << "\n";
        test_solve(slope, intercept, n_min, n_max);
    }


    return 0;
}
