#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cmath>

#define EXACT 0

const double delta = 0.5;
const double divisor = 1 / (1024.0 * 2 * delta * delta);

#if EXACT
double my_exp(double x) {
    return std::exp(-x/(2*delta*delta));
}
#else
double my_exp(double x) {
    x = 1 - x * divisor;
    x *= x; x *= x; x *= x; x *= x;
    x *= x; x *= x; x *= x; x *= x;
    x *= x; x *= x;
    return x;
}
#endif

struct Event {
  double s12;
  double s13;
  double s24;
  double s34;
  double s134;
};

const std::vector<Event> read_file(const std::string filename, const size_t n_events) {
  std::fstream file(filename, std::ios_base::in);
  if (!file) {
    std::cerr << "Error opening file " << filename << std::endl;
    exit(-1);
  }

  std::vector<Event> events;
  events.reserve(std::min((size_t) 5000000, n_events));

  std::string line;
  while (std::getline(file, line) && events.size() < n_events) {
    Event e;
    file >> e.s12 >> e.s13 >> e.s24 >> e.s34 >> e.s134;
    events.push_back(e);
  }
  return events;
}

double compute_distance(const std::vector<Event> &data_1, const std::vector<Event> &data_2, const bool upper_only) {
    double total = 0;
    #pragma omp parallel for reduction(+:total)
    for (size_t i=0; i < data_1.size(); ++i) {
        auto event_1 = data_1[i];
        for (size_t j=(upper_only ? i+1 : 0); j < data_2.size(); ++j) {
            auto event_2 = data_2[j];
            double distance_squared =
                (event_1.s12 - event_2.s12) * (event_1.s12 - event_2.s12) +
                (event_1.s13 - event_2.s13) * (event_1.s13 - event_2.s13) +
                (event_1.s24 - event_2.s24) * (event_1.s24 - event_2.s24) +
                (event_1.s34 - event_2.s34) * (event_1.s34 - event_2.s34) +
                (event_1.s134 - event_2.s134) * (event_1.s134 - event_2.s134);
            total += my_exp(distance_squared);
        }
    }
    return total;
}

int main(int argc, char *argv[]) {
  #if EXACT
    std::cout << "Running in exact mode" << std::endl;
  #else
    std::cout << "Running in approximate mode" << std::endl;
  #endif
  size_t n_events = std::numeric_limits<size_t>::max();
  if (argc == 4) {
    n_events = std::stoi(argv[3]);
    std::cout << "Limiting number of events to: " << n_events << std::endl;
  } else if (argc != 3) {
    std::cerr << "Bad arguments given, usage:" << std::endl;
    std::cerr << "    " << argv[0] << " path/to/d0_data.csv path/to/d0bar_data.csv" << std::endl;
  }

  const auto data_1 = read_file(argv[1], n_events);
  std::cout << "data_1.size() = " << data_1.size() << std::endl;

  const auto data_2 = read_file(argv[2], n_events);
  std::cout << "data_2.size() = " << data_2.size() << std::endl;

  double dist_11 = compute_distance(data_1, data_1, true);
  dist_11 /= data_1.size()*(data_1.size()-1);
  std::cout << "dist_11 = " << dist_11 << std::endl;

  double dist_22 = compute_distance(data_2, data_2, true);
  dist_22 /= data_2.size()*(data_2.size()-1);
  std::cout << "dist_22 = " << dist_22 << std::endl;

  double dist_12 = compute_distance(data_1, data_2, false);
  dist_12 /= data_1.size()*data_2.size();
  std::cout << "dist_12 = " << dist_12 << std::endl;

  std::cout << "T = " << dist_11 + dist_22 - dist_12 << std::endl;

  return 0;
}
