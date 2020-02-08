//
// Created by Logan Morrison on 2/4/20.
//

//
// Created by Logan Morrison on 2019-05-24.
//

#include "darksun/standard_model.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <ctime>
#include <cmath>

using namespace darksun;

const double log_T_min = _sm_log_temperature_min;
const double log_T_max = _sm_log_temperature_max;
const size_t num_Ts = _data_sm_sqrt_gstars.size();


class StandardModelTest : public ::testing::Test {
protected:
    void SetUp() override {
        sm_Ts.reserve(num_Ts);
        double step_log_T = (log_T_max - log_T_min) / double(num_Ts - 1);
        for (size_t n = 0; n < num_Ts; n++) {
            double log_T = log_T_min + n * step_log_T;
            sm_Ts.push_back(pow(10.0, log_T));
        }
    }

    // void TearDown() override {}
    std::vector<double> sm_Ts;
};

TEST_F(StandardModelTest, TestSqrtGStar) {
    for (size_t i = 0; i < num_Ts; i++) {
        double approx = sm_sqrt_gstar(sm_Ts[i]);
        double exact = _data_sm_sqrt_gstars[i];
        double frac_diff = fabs(approx - exact) / exact;
        ASSERT_LE(frac_diff, 5e-2);
    }
}

TEST_F(StandardModelTest, TestGeff) {
    for (size_t i = 0; i < num_Ts; i++) {
        double approx = sm_dof_energy(sm_Ts[i]);
        double exact = _data_sm_geffs[i];
        double frac_diff = fabs(approx - exact) / exact;
        ASSERT_LE(frac_diff, 5e-2);
    }
}

TEST_F(StandardModelTest, TestHeff) {
    for (size_t i = 0; i < num_Ts; i++) {
        double approx = sm_dof_entropy(sm_Ts[i]);
        double exact = _data_sm_heffs[i];
        double frac_diff = fabs(approx - exact) / exact;
        ASSERT_LE(frac_diff, 5e-2);
    }
}


int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
