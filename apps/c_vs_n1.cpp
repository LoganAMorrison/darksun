/*
 * File for generating the data for lec1=0.1 and lec2=1.0 as well as a value
 * of c such that the delta obtains the correct relic density at N=7.
 */

#include <boost/timer/progress_display.hpp>
#include <darksun/darksun.hpp>
#include <darksun/scanner.hpp>
#include <filesystem>
#include <mutex>
#include <thread>

using namespace darksun;

static constexpr size_t NUM_N = 150;
static constexpr double N_MIN = 5;
static constexpr double N_MAX = 35;
static constexpr double N_STP = (N_MAX - N_MIN) / double(NUM_N - 1);

static constexpr size_t NUM_CS = 151;
static constexpr double LOG_C_MIN = -3.0;
static constexpr double LOG_C_MAX = 1.0;
static constexpr double LOG_C_STP =
    (LOG_C_MAX - LOG_C_MIN) / double(NUM_CS - 1);

static constexpr double LEC1 = 0.1;
static constexpr double LEC2 = 1.0;
static constexpr double LAM = 1e-6;
static constexpr double XI_INF = 1e-2;

const std::string FNAME =
    "/home/logan/Research/DarkSun/cpp/rundata/c_vs_n_lam1.csv";

static boost::timer::progress_display progress(NUM_N *NUM_CS);
static std::mutex progress_mutex;

bool set_model(size_t i, DarkSunParameters &params) {
  if (i < NUM_CS * NUM_N) {
    int idx_n = i % NUM_N;
    int idx_c = (i - idx_n) / NUM_N;

    params.n = N_MIN + idx_n * N_STP;
    params.c = pow(10.0, LOG_C_MIN + idx_c * LOG_C_STP);
    params.lec1 = LEC1;
    params.lec2 = LEC2;
    params.xi_inf = XI_INF;
    params.lam = LAM;

    {
      std::lock_guard<std::mutex> lock(progress_mutex);
      ++progress;
    }

    return false;
  } else {
    return true;
  }
}

int main() {
  Scanner s(FNAME, set_model);
  s.scan();
}
