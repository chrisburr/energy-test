#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cmath>

#include <args.hxx>

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

double compute_statistic(const std::vector<Event> &data_1, const std::vector<Event> &data_2, const bool debug = false) {
  double dist_11 = compute_distance(data_1, data_1, true);
  dist_11 /= data_1.size()*(data_1.size()-1);
  if (debug)
    std::cout << "    dist_11 = " << dist_11 << std::endl;

  double dist_22 = compute_distance(data_2, data_2, true);
  dist_22 /= data_2.size()*(data_2.size()-1);
  if (debug)
    std::cout << "    dist_22 = " << dist_22 << std::endl;

  double dist_12 = compute_distance(data_1, data_2, false);
  dist_12 /= data_1.size()*data_2.size();
  if (debug)
    std::cout << "    dist_12 = " << dist_12 << std::endl;

  return dist_11 + dist_22 - dist_12;
}

int main(int argc, char *argv[]) {
  #if EXACT
    std::cout << "Running in exact mode" << std::endl;
  #else
    std::cout << "Running in approximate mode" << std::endl;
  #endif

  // Parse the command line arguments
  args::ArgumentParser parser("CPU based energy test");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::ValueFlag<size_t> n_permutations(parser, "n_permutations", "Number of permutations to run", {"n-permutations"});
  args::ValueFlag<size_t> max_events_1(parser, "max events 1", "Maximum number of events to use from dataset 1", {"max-events-1"});
  args::ValueFlag<size_t> max_events_2(parser, "max events 2", "Maximum number of events to use from dataset 2", {"max-events-2"});
  args::ValueFlag<size_t> max_events(parser, "max events", "Max number of events in each dataset", {"max-events"});
  args::Positional<std::string> filename_1(parser, "dataset 1", "Filename for the first dataset");
  args::Positional<std::string> filename_2(parser, "dataset 2", "Filename for the second dataset");
  try {
    parser.ParseCLI(argc, argv);
    if (!filename_1 || !filename_2)
      throw args::ParseError("Two dataset filenames must be given");
    if ((max_events_1 || max_events_2) && max_events)
      throw args::ParseError("--max-events cannot be used with --max-events-1 or --max-events-2");
  } catch (args::Help) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  // Parse the maximum number of events to use
  size_t data_1_limit = std::numeric_limits<size_t>::max();
  size_t data_2_limit = std::numeric_limits<size_t>::max();
  if (max_events) {
    data_1_limit = args::get(max_events);
    data_2_limit = args::get(max_events);
  } else {
    if (max_events_1)
      data_1_limit = args::get(max_events_1);
    if (max_events_2)
      data_2_limit = args::get(max_events_2);
  }

  // Load the data
  const auto dataset_1 = read_file(args::get(filename_1), data_1_limit);
  std::cout << "Dataset 1 size is " << dataset_1.size() << std::endl;
  const auto dataset_2 = read_file(args::get(filename_2), data_2_limit);
  std::cout << "Dataset 2 size is " << dataset_2.size() << std::endl;

  std::cout << std::endl;

  // Compute the test statistic for the current dataset
  std::cout << "Calculating test statistic for nominal dataset:" << std::endl;
  const double real_test_statistic = compute_statistic(dataset_1, dataset_2, true);
  std::cout << "  T = " << real_test_statistic << std::endl;

  std::cout << std::endl;

  if (n_permutations) {
    // Merge the vectors of events so we can shuffle them
    std::vector<Event> all_events;
    all_events.insert(all_events.end(), dataset_1.begin(), dataset_1.end());
    all_events.insert(all_events.end(), dataset_2.begin(), dataset_2.end());

    auto N = args::get(n_permutations);
    for (size_t i = 0; i < N; ++i) {
      std::cout << "Calculating permutation " << i+1 << " of " << N;
      std::random_shuffle(all_events.begin(), all_events.end());

      const std::vector<Event> data_1(all_events.begin(), all_events.begin()+dataset_1.size());
      const std::vector<Event> data_2(all_events.begin()+dataset_1.size(), all_events.end());

      const double test_statistic = compute_statistic(data_1, data_2);
      std::cout << ": T=" << test_statistic << std::endl;
    }
  }
  return 0;
}
