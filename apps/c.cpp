
#include <boost/timer/progress_display.hpp>
#include <darksun/darksun.hpp>
#include <filesystem>
#include <fmt/core.h>
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

static constexpr std::array<double, 4> CS = {0.01, 0.1, 0.5, 1.5};
static const std::array<std::string, CS.size()> FNAMES = {
    std::filesystem::current_path().append("../rundata/c=0.01.csv"),
    std::filesystem::current_path().append("../rundata/c=0.1.csv"),
    std::filesystem::current_path().append("../rundata/c=0.5.csv"),
    std::filesystem::current_path().append("../rundata/c=1.5.csv")};

static boost::timer::progress_display progress(CS.size() * NUM_N * NUM_LAM);
static std::mutex progress_mutex;

bool set_model(size_t i, DarkSunParameters &params, double c) {
  if (i < NUM_LAM * NUM_N) {
    int idx_n = i % NUM_N;
    int idx_lam = (i - idx_n) / NUM_N;

    params.n = N_MIN + idx_n * N_STP;
    params.lam = pow(10.0, LOG_LAM_MIN + idx_lam * LOG_LAM_STP);
    params.c = c;

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

  std::array<std::thread, CS.size()> threads;

  for (size_t i = 0; i < CS.size(); i++) {
    auto f = [i](int counter, DarkSunParameters &params) {
      return set_model(counter, params, CS[i]);
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
