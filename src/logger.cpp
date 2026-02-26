#include "logger.h"

#include <mpi.h>
#include <spdlog/common.h>

#include <csignal>
#include <iostream>
#include <vector>

#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"

namespace dendro {
namespace logger {

std::shared_ptr<spdlog::logger> dendrologger;

void crash_signal_log_flush_handler(int signum) {
    fprintf(stderr, "Caught signal %d, Flushing logs before exit!", signum);
    if (dendro::logger::dendrologger) {
        dendrologger->flush();
    }
    signal(signum, SIG_DFL);
    raise(signum);
}

void setup_crash_handler() {
    // segmentation fault handling
    signal(SIGSEGV, crash_signal_log_flush_handler);
    // abort signal handling
    signal(SIGABRT, crash_signal_log_flush_handler);
    // illegal instruction handling
    signal(SIGILL, crash_signal_log_flush_handler);
    // floating point exception handling
    signal(SIGFPE, crash_signal_log_flush_handler);

    // termination signals
    // sig interrupt
    signal(SIGINT, crash_signal_log_flush_handler);
    // sig
    signal(SIGTERM, crash_signal_log_flush_handler);
}

void initialize(const std::string& logger_name, const std::string& log_filename,
                int logging_rank, spdlog::level::level_enum file_level,
                spdlog::level::level_enum console_level,
                bool force_debug_flush) {
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank == logging_rank) {
        std::cout << "NOW BUILDING LOGGER ON RANK " << logging_rank << " - "
                  << world_rank << std::endl;
        try {
            std::vector<spdlog::sink_ptr> sinks;

            // create the console sink
            auto console_sink =
                std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
            console_sink->set_level(console_level);
            sinks.push_back(console_sink);

            // create the file sink
            auto file_sink =
                std::make_shared<spdlog::sinks::basic_file_sink_mt>(
                    log_filename, true);
            file_sink->set_level(file_level);
            sinks.push_back(file_sink);

            // create the logger
            dendrologger = std::make_shared<spdlog::logger>(
                logger_name, begin(sinks), end(sinks));

            // set the overall minimum level to the most verbose setting
            dendrologger->set_level(spdlog::level::trace);

            // make sure warn or higher forces a flush
            if (force_debug_flush) {
                dendrologger->flush_on(spdlog::level::debug);
            } else {
                dendrologger->flush_on(spdlog::level::warn);
            }

            // make this logger the default
            spdlog::set_default_logger(dendrologger);
        } catch (const spdlog::spdlog_ex& ex) {
            std::cerr << "Log initialization failed: " << ex.what()
                      << std::endl;
        }
    }
}

void initialize(const std::string& logger_name, const std::string& log_filename,
                int logging_rank, int file_level, int console_level,
                bool force_debug_flush) {
    initialize(logger_name, log_filename, logging_rank,
               static_cast<spdlog::level::level_enum>(file_level),
               static_cast<spdlog::level::level_enum>(console_level),
               force_debug_flush);
}

}  // namespace logger
}  // namespace dendro
