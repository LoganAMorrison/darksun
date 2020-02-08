//
// Created by Logan Morrison on 2/4/20.
//

#include "darksun/models.hpp"
#include <gtest/gtest.h>

using namespace darksun;

TEST(Test1, Test11) {
    double exact1 = 7.589746948024013e-18;
    double exact2 = 7.07413e-40;
    double exact3 = 3.40632e-67;

    DarkSUNModel model1{1.0, 5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, false};
    DarkSUNModel model2{1.0, 10, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, false};
    DarkSUNModel model3{1.0, 15, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, false};

    double approx1 = model1.thermal_cross_section_2eta_2delta(1.0);
    double approx2 = model2.thermal_cross_section_2eta_2delta(1.0);
    double approx3 = model3.thermal_cross_section_2eta_2delta(1.0);

    double frac_diff1 = fabs(approx1 - exact1) / exact1;
    double frac_diff2 = fabs(approx2 - exact2) / exact2;
    double frac_diff3 = fabs(approx3 - exact3) / exact3;

    ASSERT_LE(frac_diff1, 1e-3);
    ASSERT_LE(frac_diff2, 1e-3);
    ASSERT_LE(frac_diff3, 1e-3);
}

TEST(Test1, Test12) {
    DarkSUNModel model1{1.0, 10, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, false};
    DarkSUNModel model2{1.0, 20, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, false};
    DarkSUNModel model3{1.0, 30, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, false};

    double approx1 = model1.thermal_cross_section_2eta_4eta(1.0);
    double approx2 = model2.thermal_cross_section_2eta_4eta(1.0);
    double approx3 = model3.thermal_cross_section_2eta_4eta(1.0);

    std::cout << "approx1 = " << approx1 << std::endl;
    std::cout << "approx2 = " << approx2 << std::endl;
    std::cout << "approx3 = " << approx3 << std::endl;
}


int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
