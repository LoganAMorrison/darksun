
#include <boost/timer/progress_display.hpp>
#include <darksun/darksun.hpp>
#include <darksun/scanner.hpp>
#include <filesystem>
#include <mutex>
#include <thread>

using namespace darksun;

static constexpr size_t N_MIN = 5;
static constexpr size_t N_MAX = 100;
static constexpr size_t NUM_N = N_MAX - N_MIN + 1;
static constexpr size_t N_STP = 1;

static constexpr size_t NUM_LAM = 100;
static constexpr double LOG_LAM_MIN = -7.0;
static constexpr double LOG_LAM_MAX = 1.0;
static constexpr double LOG_LAM_STP =
    (LOG_LAM_MAX - LOG_LAM_MIN) / double(NUM_LAM - 1);

const std::string FNAME =
    std::filesystem::current_path().append("../rundata/bm.csv");

static boost::timer::progress_display progress(NUM_N *NUM_LAM);
static std::mutex progress_mutex;

bool set_model(size_t i, DarkSunParameters &params) {
  if (i < NUM_LAM * NUM_N) {
    int idx_n = i % NUM_N;
    int idx_lam = (i - idx_n) / NUM_N;

    params.n = N_MIN + idx_n * N_STP;
    params.lam = pow(10.0, LOG_LAM_MIN + idx_lam * LOG_LAM_STP);

    ++progress;

    return false;
  } else {
    return true;
  }
}

int main() {
  Scanner s(FNAME, set_model);
  s.scan();
}
