#include <boost/timer/progress_display.hpp>
#include <darksun/darksun.hpp>
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

static constexpr std::array<double, 2> LEC1S = {0.01, 1.0};
static const std::array<std::string, LEC1S.size()> FNAMES = {
    std::filesystem::current_path().append("../rundata/lec1=0.01.csv"),
    std::filesystem::current_path().append("../rundata/lec1=1.0.csv")};

static boost::timer::progress_display progress(LEC1S.size() * NUM_N * NUM_LAM);
static std::mutex progress_mutex;

bool set_model(size_t i, DarkSunParameters &params, double lec1) {
  if (i < NUM_LAM * NUM_N) {
    int idx_n = i % NUM_N;
    int idx_lam = (i - idx_n) / NUM_N;

    params.n = N_MIN + idx_n * N_STP;
    params.lam = pow(10.0, LOG_LAM_MIN + idx_lam * LOG_LAM_STP);
    params.lec1 = lec1;

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

  std::array<std::thread, LEC1S.size()> threads;

  for (size_t i = 0; i < LEC1S.size(); i++) {
    auto f = [i](int counter, DarkSunParameters &params) {
      return set_model(counter, params, LEC1S[i]);
    };
    Scanner s(FNAMES[i], f);
    threads[i] = std::thread(&Scanner::scan, s);
  }
  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }
}
