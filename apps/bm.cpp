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

static constexpr size_t NUM_LAM = 150;
static constexpr double LOG_LAM_MIN = -7.0;
static constexpr double LOG_LAM_MAX = 1.0;
static constexpr double LOG_LAM_STP =
    (LOG_LAM_MAX - LOG_LAM_MIN) / double(NUM_LAM - 1);

static constexpr double LEC1 = 0.1;
static constexpr double LEC2 = 1.0;

// Value of c that yields correct delta RD at N = 7
static constexpr double C = 0.666544284531189;

const std::string FNAME = std::filesystem::current_path().append(
    "../rundata/bm_lec1=0.1_lec2=1.0.csv");

static boost::timer::progress_display progress(NUM_N *NUM_LAM);
static std::mutex progress_mutex;

bool set_model(size_t i, DarkSunParameters &params) {
  if (i < NUM_LAM * NUM_N) {
    int idx_n = i % NUM_N;
    int idx_lam = (i - idx_n) / NUM_N;

    params.n = N_MIN + idx_n * N_STP;
    params.lam = pow(10.0, LOG_LAM_MIN + idx_lam * LOG_LAM_STP);
    params.lec1 = LEC1;
    params.lec2 = LEC2;
    params.c = C;
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
