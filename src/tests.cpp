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


void test_AMCDF_HaloICA1() {
    printf("\n=== AMCDF TESTS (HALO ICA 1) ===\n");

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
    printf("\n--- Testing Boundary Behavior  ---\n");
    double extreme_low = -1000.0;
    double val_low = t.eval(extreme_low);
    double extreme_high = 1000.0;
    double val_high = t.eval(extreme_high);
    
    printf("  x=%.1f -> %.6f,   x=%.1f -> %.3e\n", extreme_low, val_low, extreme_high, val_high);
    
    TEST_CASE(std::isfinite(val_low) && val_low > 1e-2, val_low, "CDF should handle extrapolation for low values");
    TEST_CASE(std::isfinite(val_high) && val_high < 1e-8, val_high, "CDF should handle extrapolation for high values");
    TEST_CASE(val_low >= 0.0 && val_high >= 0.0, val_low, "Extrapolated densities should remain non-negative");

}

void test_AMCDF_GalaxyAbsMag() {

    // Now repeat for magnitudes file
    AMCDF t(GAL_MAGR_DENSITY);
    std::vector<double> x_vals(20);

    printf("\n--- Testing Spline Smoothness and Monotonicity (Galaxy abs mag) ---\n");
    for (int i = 0; i < 20; i++)
        x_vals[i] = -23.0 + i * 0.1; // Sample from -23 to -21
    double prev_val = -9999;
    bool is_smooth = true;
    for (double x : x_vals) {
        double val = t.eval(x);            
        TEST_CASE(std::isfinite(val), val, "All interpolated values should be finite");

        if (prev_val > -9999) {
            // Monotonic Test
            TEST_CASE(val >= prev_val, x, "Cumulative density should get larger as we increase x"); // Summed from small values up for magnitudes since negative is brighter

            // Check that moving by 0.1 in this space doesn't cause huge jumps in density (which would indicate a bad spline fit)
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
    printf("\n--- Testing Boundary Behavior (Galaxy abs mag) ---\n");
    double extreme_low = -30.0;
    double val_low = t.eval(extreme_low);
    double extreme_high = -5.0;
    double val_high = t.eval(extreme_high);

    TEST_CASE(std::isfinite(val_low) && val_low < 1e-8, val_low, "CDF should handle extrapolation for low values");
    TEST_CASE(std::isfinite(val_high) && val_high > 1e-2, val_high, "CDF should handle extrapolation for high values");
    TEST_CASE(val_low >= 0.0 && val_high >= 0.0, val_low, "Extrapolated densities should remain non-negative");



    printf("\n--- Testing Inverted Spline Smoothness and Monotonicity (Galaxy abs mag) ---\n");
    std::vector<double> densities(100);
    for (int i = 0; i < 100; i++) {
        densities[i] = 1e-8 * pow(10, i * 0.1); // Sample densities logarithmically
    }
    prev_val = -9999;
    is_smooth = true;
    bool all_finite = true;
    for (double d : densities) {
        double mag = t.eval_inverse(d);       
        if (!std::isfinite(mag)) {
            all_finite = false;
            printf("    WARNING: Non-finite magnitude for density %.2e\n", d);
        }     

        if (prev_val > -9999 && all_finite) {
            // Monotonic Test
            TEST_CASE(mag >= prev_val, mag, "Mag should get dimmer (larger number) as density goes up"); 

            double ratio = mag / prev_val;
            if (ratio > 1.25 || ratio < 0.75) {
                is_smooth = false;
                printf("    WARNING: Large jump detected (ratio=%.2e)\n", ratio);
            }
        }
        prev_val = mag;
    }
    TEST_CASE(all_finite, all_finite, "All magnitudes from inverse spline should be finite");
    TEST_CASE(is_smooth, is_smooth, "Spline should be reasonably smooth (no huge jumps)");
}


void test_AMCDF_GalaxyGmR() {
    std::cout << "\n=== AMCDF TESTS (GALAXY COLOR g-r) ===\n";

    AMCDF t(GAL_COLOR_GMR_DENSITY);
    std::vector<double> x_vals(20);

    printf("\n--- Testing Spline Smoothness and Monotonicity (Galaxy g-r) ---\n");
    for (int i = 0; i < 20; i++)
        x_vals[i] = 0.2 + i * 0.04;
    double prev_val = -9999;
    bool is_smooth = true;
    for (double x : x_vals) {
        double val = t.eval(x);            
        TEST_CASE(std::isfinite(val), val, "All interpolated values should be finite");

        if (prev_val > -9999) {
            // Monotonic Test
            TEST_CASE(val <= prev_val, x, "Cumulative density should get smaller as we increase x"); // Summed from high g-r (redder) down

            // Check that moving by 0.1 in this space doesn't cause huge jumps in density (which would indicate a bad spline fit)
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
    printf("\n--- Testing Boundary Behavior (Galaxy g-r) ---\n");
    double extreme_low = -5.0;
    double val_low = t.eval(extreme_low);
    double extreme_high = 5.0;
    double val_high = t.eval(extreme_high);

    TEST_CASE(std::isfinite(val_low) && val_low > 1e-2, val_low, "CDF should handle extrapolation for low values");
    TEST_CASE(std::isfinite(val_high) && val_high < 1e-8, val_high, "CDF should handle extrapolation for high values");
    TEST_CASE(val_low >= 0.0 && val_high >= 0.0, val_low, "Extrapolated densities should remain non-negative");

}

void test_GalaxyMagMatcher() {
    printf("\n=== GalaxyMagMatcher TESTS ===\n");

    auto matcher = GalaxyMagMatcher::get();

    // The sim box will start at 1.56e-8 so we need at least this low density.
    // Currently cannot handle much higher densities than 0.05
    std::vector<double> test_densities = {1e-9, 1e-8, 1e-6, 1e-4, 0.01, 0.03, 0.05}; 
    for (double dens : test_densities) {
        double mag = matcher.match(dens);
        TEST_CASE(std::isfinite(mag) && mag < -5 && mag > -30 , mag, "Matched magnitude should be finite and reasonable");
        printf("  Density=%.2e -> Mag=%.3f\n", dens, mag);
    }

    // Test that higher density gives fainter magnitude (higher number in mags)
    double mag1 = matcher.match(1e-4);
    double mag2 = matcher.match(1e-3);
    TEST_CASE(mag2 > mag1, mag2 - mag1, "Higher density should give fainter magnitude");
}

int main () {
    test_TabulatedDensityFunction();
    test_AMCDF_HaloICA1();
    test_AMCDF_GalaxyAbsMag();
    test_AMCDF_GalaxyGmR();
    test_GalaxyMagMatcher();

    std::cout << "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" << std::endl;
    std::cout << "All tests passed successfully." << std::endl;
    std::cout << "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" << std::endl;
    return 0;
}