

#include <filesystem>
#include <iostream>
#include <thread>

#include <framework/utils/argparse.h>

#include "Task.h"

using namespace std::literals;

std::ostream& printUsage(std::ostream& out) {
  return out
         << R"""(usage: task1 [-b <max-bounces>] [-t [<num-threads>]|--num-threads <num-threads>] <config-file>
	-w --res-x              override horizontal resolution
	-h --res-y              override vertical resolution
	-b                      set maximum number of bounces
	-t --num-threads        render in parallel using the specified number of threads (default when only -t is given: )"""
         << std::thread::hardware_concurrency() << R"""()
	-h --help               display this message
)""";
}

int taskMain(int argc, char* argv[]) {
  int res_x = 0;
  int res_y = 0;
  int num_threads = 1;
  std::filesystem::path input_file;
  bool no_cache = true;
  int max_bounces = 3;

  for (const char* const* a = argv + 1; *a; ++a) {
    if (parseIntArgument(num_threads, a, "t"sv,
                         static_cast<int>(std::thread::hardware_concurrency())))
      ;
    else if (parseIntArgument(num_threads, a, "-num-threads"sv))
      ;
    else if (parseIntArgument(res_x, a, "w"sv) ||
             parseIntArgument(res_x, a, "-res-x"sv))
      ;
    else if (parseIntArgument(res_y, a, "h"sv) ||
             parseIntArgument(res_y, a, "-res-y"sv))
      ;
    else if (parseIntArgument(max_bounces, a, "b"sv))
      ;
    else if (parseBoolFlag(a, "h"sv) || parseBoolFlag(a, "-help"sv))
      printUsage(std::cout) << std::endl;
    else if (input_file.empty())
      input_file = *a;
    else
      throw usage_error("invalid argument");
  }

  if (input_file.empty()) throw usage_error("missing config file argument");

  std::filesystem::path output_dir = "output/task1";

  std::filesystem::create_directories(output_dir);

  Task task(output_dir, input_file, no_cache);
  const std::filesystem::__cxx11::path fn = (output_dir / input_file.filename()).replace_extension(".png");
  task.render(fn,
            res_x, res_y, max_bounces, num_threads);

  return 0;
}
