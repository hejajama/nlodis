#include "unit_test_framework.hpp"
#include "../dipole/bkdipole/interpolation.hpp"
#include <cmath>

TEST(interpolate_square)
{
    double x[5] = {0, 1, 2, 3, 4};
    double y[5] = {0, 1, 4, 9, 16};
    Interpolator inter(x,y,5);
    inter.SetMethod(InterpolationMethod::SPLINE);
    inter.Initialize();
    
    ASSERT_ALMOST_EQUAL(inter.Evaluate(0),0,1e-6);
    ASSERT_ALMOST_EQUAL(inter.Evaluate(1),1,1e-6);
    ASSERT_ALMOST_EQUAL(inter.Evaluate(2),4,1e-6);
    ASSERT_ALMOST_EQUAL(inter.Evaluate(3),9,1e-6);
    ASSERT_ALMOST_EQUAL(inter.Evaluate(4),16,1e-6);

    // Same test with std::vector
    std::vector<double> xvec(x, x + 5);
    std::vector<double> yvec(y, y + 5);
    Interpolator inter2(xvec, yvec);
    inter2.SetMethod(InterpolationMethod::SPLINE);
    inter2.Initialize();
    ASSERT_ALMOST_EQUAL(inter2.Evaluate(0),0,1e-6);
    ASSERT_ALMOST_EQUAL(inter2.Evaluate(1),1,1e-6);
    ASSERT_ALMOST_EQUAL(inter2.Evaluate(2),4,1e-6);
    ASSERT_ALMOST_EQUAL(inter2.Evaluate(3),9,1e-6);
    ASSERT_ALMOST_EQUAL(inter2.Evaluate(4),16,1e-6);
}

TEST(interpolate_out_of_range)
{
    double x[5] = {0, 1, 2, 3, 4};
    double y[5] = {0, 1, 4, 9, 16};
    Interpolator inter(x,y,5);
    inter.SetMethod(InterpolationMethod::SPLINE);
    inter.SetFreeze(true); // Freeze below and above the given data
    inter.Initialize();
    
    ASSERT_ALMOST_EQUAL(inter.Evaluate(-1),0,1e-6);
    ASSERT_ALMOST_EQUAL(inter.Evaluate(5),16,1e-6);

    Interpolator inter2(x,y,5);
    inter2.SetMethod(InterpolationMethod::SPLINE);
    try {
        inter2.Evaluate(-1);
        ASSERT_TRUE(false); // Should not reach here
    } catch (const std::out_of_range& e) {
        ASSERT_TRUE(true); // Exception caught as expected
    }
    
}

TEST(interpolate_copy_constructor)
{
    std::vector<double> x = {0, 1, 2, 3, 4};
    std::vector<double> y = {0, 1, 4, 9, 16};
    Interpolator inter(x, y);
    inter.Initialize();

    Interpolator inter2 (inter);

    inter.Clear();
    inter.Clear();

    ASSERT_ALMOST_EQUAL(inter2.Evaluate(3),9,1e-6);

    Interpolator* inter3 = new Interpolator(x,y);

    Interpolator inter4(*inter3);

    inter3->Clear();
    delete inter3;

    ASSERT_ALMOST_EQUAL(inter4.Evaluate(3),9,1e-6);

    // Same test but construct interpolator using x and y arrays, as 
    // in this case Interpolator does not reserve memory for x and y data
    double x_arr[5] = {0, 1, 2, 3, 4};
    double y_arr[5] = {0, 1, 4, 9, 16};
    Interpolator inter5(x_arr, y_arr, 5);
    inter5.Initialize();

    Interpolator inter6(inter5);

    inter5.Clear();

    ASSERT_ALMOST_EQUAL(inter6.Evaluate(3), 9, 1e-6);


}

TEST(interpolate_assignment_operator)
{
    std::vector<double> x = {0, 1, 2, 3, 4};
    std::vector<double> y = {0, 1, 4, 9, 16};
    Interpolator inter(x, y);
    inter.Initialize();

    Interpolator inter2;

    inter2 = inter;


    inter.Clear();
    ASSERT_ALMOST_EQUAL(inter2.Evaluate(3), 9, 1e-6);

    Interpolator* inter3 = new Interpolator(x, y);
    inter3->Initialize();

    Interpolator inter4;
    inter4 = *inter3;

    inter3->Clear();
    delete inter3;

    ASSERT_ALMOST_EQUAL(inter4.Evaluate(3), 9, 1e-6);

    // Same test but construct interpolator using x and y arrays
    double x_arr[5] = {0, 1, 2, 3, 4};
    double y_arr[5] = {0, 1, 4, 9, 16};
    Interpolator inter5(x_arr, y_arr, 5);
    inter5.Initialize();

    Interpolator inter6;
    inter6 = inter5;

    inter5.Clear();

    ASSERT_ALMOST_EQUAL(inter6.Evaluate(3), 9, 1e-6);
}

// ---- Log-grid interpolation tests ----

TEST(interpolate_log_power_law)
{
    // f(x) = x^2 is linear in log-log space: log(f) = 2*log(x)
    // A cubic spline on the log-log grid should reproduce it exactly.
    std::vector<double> x, y;
    for (int i = 0; i < 6; ++i) {
        double xi = std::pow(10.0, i - 1); // 0.1, 1, 10, 100, 1000, 10000
        x.push_back(xi);
        y.push_back(xi * xi);
    }
    Interpolator inter(x, y, true);

    // At grid points the result must match exactly
    for (int i = 0; i < 6; ++i) {
        ASSERT_ALMOST_EQUAL(inter.Evaluate(x[i]), y[i], 1e-6);
    }

    // Between grid points the log-log spline should still be highly accurate
    ASSERT_ALMOST_EQUAL(inter.Evaluate(0.5),   0.25,   1e-6);
    ASSERT_ALMOST_EQUAL(inter.Evaluate(3.0),   9.0,    1e-4);
    ASSERT_ALMOST_EQUAL(inter.Evaluate(50.0),  2500.0, 1e-2);
}

TEST(interpolate_log_invalid_construction_x)
{
    // Negative x value must throw during construction
    {
        std::vector<double> x_neg = {-1.0, 1.0, 2.0};
        std::vector<double> y_neg = {1.0, 1.0, 4.0};
        try {
            Interpolator inter(x_neg, y_neg, true);
            ASSERT_TRUE(false); // should not reach here
        } catch (const std::invalid_argument&) {
            ASSERT_TRUE(true);
        }
    }

    // Zero x value must throw during construction
    {
        std::vector<double> x_zero = {0.0, 1.0, 2.0};
        std::vector<double> y_zero = {1.0, 1.0, 4.0};
        try {
            Interpolator inter2(x_zero, y_zero, true);
            ASSERT_TRUE(false);
        } catch (const std::invalid_argument&) {
            ASSERT_TRUE(true);
        }
    }
}

TEST(interpolate_log_invalid_construction_y)
{
    // Zero y value must throw during construction
    {
        std::vector<double> x = {1.0, 2.0, 3.0};
        std::vector<double> y_zero = {1.0, 0.0, 4.0};
        try {
            Interpolator inter(x, y_zero, true);
            ASSERT_TRUE(false);
        } catch (const std::invalid_argument&) {
            ASSERT_TRUE(true);
        }
    }

    // Negative y value must throw during construction
    {
        std::vector<double> x = {1.0, 2.0, 3.0};
        std::vector<double> y_neg = {1.0, -3.0, 4.0};
        try {
            Interpolator inter(x, y_neg, true);
            ASSERT_TRUE(false);
        } catch (const std::invalid_argument&) {
            ASSERT_TRUE(true);
        }
    }
}

TEST(interpolate_log_eval_nonpositive_x)
{
    // Evaluating at x <= 0 must throw when interpolate_log=true
    std::vector<double> x = {1.0, 2.0, 4.0, 8.0};
    std::vector<double> y = {1.0, 4.0, 16.0, 64.0};
    Interpolator inter(x, y, true);

    try {
        inter.Evaluate(-1.0);
        ASSERT_TRUE(false);
    } catch (const std::exception&) {
        ASSERT_TRUE(true);
    }

    try {
        inter.Evaluate(0.0);
        ASSERT_TRUE(false);
    } catch (const std::exception&) {
        ASSERT_TRUE(true);
    }
}

TEST(interpolate_log_derivative)
{
    // f(x) = x^2  =>  f'(x) = 2x,  f''(x) = 2
    std::vector<double> x, y;
    for (int i = 0; i < 6; ++i) {
        double xi = std::pow(10.0, i - 1);
        x.push_back(xi);
        y.push_back(xi * xi);
    }
    Interpolator inter(x, y, true);

    ASSERT_ALMOST_EQUAL(inter.Derivative(1.0),  2.0,  1e-4);
    ASSERT_ALMOST_EQUAL(inter.Derivative(10.0), 20.0, 1e-3);

    ASSERT_ALMOST_EQUAL(inter.Derivative2(1.0),  2.0, 1e-3);
    ASSERT_ALMOST_EQUAL(inter.Derivative2(10.0), 2.0, 1e-1);
}

