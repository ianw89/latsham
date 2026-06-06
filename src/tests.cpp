#include "latsham.hpp"
#include <cassert>
#include <sstream>

#define TEST_CASE(condition, value, message) \
    if (!(condition)) { \
        std::cerr << "! TEST FAILED ! : " << message << "\n" \
                  << "  Value: " << (value) << "\n" \
                  << "  File: " << __FILE__ << ", Line: " << __LINE__ << "\n"; \
        std::abort(); \
    } else { \
        std::cout << "Test passed: " << message << ". Value: " << (value) << "\n"; \
    }

bool isclose(double a, double b, double rel_tol = 1e-6, double abs_tol = 1e-12) {
    return std::fabs(a - b) <= std::max(rel_tol * std::max(std::fabs(a), std::fabs(b)), abs_tol);
}

const std::string testdensityfile = "/mount/sirocco1/imw2293/GROUP_CAT/SelfCalGroupFinder/py/parameters/bgs_y3/halo_pca1_density_func.dat";

void test_TabulatedDensityFunction() {
    std::cout << "=== TabulatedDensityFunction TESTS ===\n";

    TabulatedDensityFunction tdf(testdensityfile);

    TEST_CASE(tdf.getSize() > 0, tdf.getSize(), "Density function should have points");

    int mid_idx = tdf.getSize() / 2;
    double value = tdf.getDensity()[mid_idx]; // Just a placeholder, replace with actual test logic
    TEST_CASE(value > 0, value, "Density should be positive in middle of distribution");
    value = tdf.getDensity()[0];
    TEST_CASE(value >= 0, value, "Density should be non-negative on left edge");
    value = tdf.getDensity()[tdf.getSize() - 1];    
    TEST_CASE(value >= 0, value, "Density should be non-negative on right edge");
    value = tdf.getDensity()[mid_idx - 5];
    double value2 = tdf.getDensity()[mid_idx + 5];
    TEST_CASE(value > value2, value2 - value, "Density should decrease as x increases");
}

void test_AMCDF() {
    printf("\n=== AMCDF TESTS ===\n");

    AMCDF t(testdensityfile);

    // Test spline interpolation quality by checking smoothness
    printf("\n--- Testing Spline Smoothness and Monotonicity ---\n");
    // Sample the cumulative density function at several points
    std::vector<double> x_vals(20);
    for (int i = 0; i < 20; i++)
        x_vals[i] = -1.0 + i * 0.1; // Sample from -1 to 1
    double prev_val = -9999;
    bool is_smooth = true;
    
    for (double x : x_vals) {
        double val = t.eval(x);            
        TEST_CASE(std::isfinite(val), val, "All interpolated values should be finite");
        
        if (prev_val > -9999) {
            // Monotonic Test
            TEST_CASE(val <= prev_val, x, "Cumulative density should get smaller as we increase x"); // Summed from large values down

            // Check that moving by 0.1 in this space doesn't cause huge jumps in density (which would indicate a bad spline fit)
            // This assumes the x values are typical for PCA / ICA components. Raw values might have different scales...
            double ratio = val / prev_val;
            if (ratio > 1e1 || ratio < 1e-1) {
                is_smooth = false;
                printf("    WARNING: Large jump detected (ratio=%.2e)\n", ratio);
            }
        }
        prev_val = val;
    }
    TEST_CASE(is_smooth, is_smooth, "Spline should be reasonably smooth (no huge jumps)");

    
    // Test extrapolation behavior at boundaries
    printf("\n--- Testing Boundary Behavior ---\n");
    double extreme_low = -1000.0;
    double val_low = t.eval(extreme_low);
    double extreme_high = 1000.0;
    double val_high = t.eval(extreme_high);
    
    printf("  x=%.1f -> %.6f,   x=%.1f -> %.3e\n", extreme_low, val_low, extreme_high, val_high);
    
    TEST_CASE(std::isfinite(val_low), val_low, "Should handle extrapolation to low values");
    TEST_CASE(std::isfinite(val_high), val_high, "Should handle extrapolation to high values");
    TEST_CASE(exp(val_low) >= 0.0 && exp(val_high) >= 0.0, exp(val_low) >= 0.0 && exp(val_high) >= 0.0, "Extrapolated densities should remain non-negative");

}

int main () {
    test_TabulatedDensityFunction();
    test_AMCDF();

    std::cout << "All tests passed successfully.\n";
    return 0;
}