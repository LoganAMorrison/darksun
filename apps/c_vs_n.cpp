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

static constexpr size_t NUM_CS = 100;
static constexpr double LOG_C_MIN = -3.0;
static constexpr double LOG_C_MAX = 1.0;
static constexpr double LOG_C_STP =
    (LOG_C_MAX - LOG_C_MIN) / double(NUM_CS - 1);

static constexpr std::array<double, 5> LAMS = {1e-7, 1e-5, 1e-3, 1e-1, 10.0};
static const std::array<std::string, LAMS.size()> FNAMES = {
    std::filesystem::current_path().append("../rundata/c_vs_n_lam=1e-7.csv"),
    std::filesystem::current_path().append("../rundata/c_vs_n_lam=1e-5.csv"),
    std::filesystem::current_path().append("../rundata/c_vs_n_lam=1e-3.csv"),
    std::filesystem::current_path().append("../rundata/c_vs_n_lam=1e-1.csv"),
    std::filesystem::current_path().append("../rundata/c_vs_n_lam=10.0.csv")};

static boost::timer::progress_display progress(LAMS.size() * NUM_N * NUM_CS);
static std::mutex progress_mutex;

bool set_model(size_t i, DarkSunParameters &params, double lam) {
  if (i < NUM_CS * NUM_N) {
    int idx_n = i % NUM_N;
    int idx_c = (i - idx_n) / NUM_N;

    params.n = N_MIN + idx_n * N_STP;
    params.c = pow(10.0, LOG_C_MIN + idx_c * LOG_C_STP);
    params.lam = lam;

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

  std::array<std::thread, LAMS.size()> threads;

  for (size_t i = 0; i < LAMS.size(); i++) {
    auto f = [i](int counter, DarkSunParameters &params) {
      return set_model(counter, params, LAMS[i]);
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
