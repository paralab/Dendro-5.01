#pragma once

#include <memory>
#include <string>

#include "spdlog/spdlog.h"

namespace dendro {
namespace logger {

// pointer to the logger that functions will use and call
extern std::shared_ptr<spdlog::logger> dendrologger;

void setup_crash_handler();

// NOTE: by calling this function we cannot capture the signals with other
// functions!
void crash_signal_log_flush_handler(int signum);

/**
 * @brief Initializes the global logger (for specific MPI rank)
 *
 * @param logger_name Name of the logger
 * @param log_filename The file to log to
 * @param logging_rank MPI Rank that performs the logging
 * @param file_level Minimum level to log to the file
 * @param console_level Minimum level to log to the console
 */
void initialize(const std::string& logger_name, const std::string& log_filename,
                int logging_rank,
                spdlog::level::level_enum file_level    = spdlog::level::debug,
                spdlog::level::level_enum console_level = spdlog::level::info,
                bool force_debug_flush                  = false);

/**
 * @brief Initializes the global logger (for specific MPI rank)
 *
 * @param logger_name Name of the logger
 * @param log_filename The file to log to
 * @param logging_rank MPI Rank that performs the logging
 * @param file_level Minimum level to log to the file
 * @param console_level Minimum level to log to the console
 *
 * @note 1 is DEBUG and 2 is INFO
 */
void initialize(const std::string& logger_name, const std::string& log_filename,
                int logging_rank, int file_level = 1, int console_level = 2,
                bool force_debug_flush = false);

struct Scope {
    std::string_view name;
};

template <typename... Args>
void trace(spdlog::format_string_t<Args...> fmt, Args&&... args) {
    if (dendrologger) {
        dendrologger->trace(fmt, std::forward<Args>(args)...);
    }
}

template <typename... Args>
void debug(spdlog::format_string_t<Args...> fmt, Args&&... args) {
    if (dendrologger) {
        dendrologger->debug(fmt, std::forward<Args>(args)...);
    }
}

template <typename... Args>
void info(spdlog::format_string_t<Args...> fmt, Args&&... args) {
    if (dendrologger) {
        dendrologger->info(fmt, std::forward<Args>(args)...);
    }
}

template <typename... Args>
void warn(spdlog::format_string_t<Args...> fmt, Args&&... args) {
    if (dendrologger) {
        dendrologger->warn(fmt, std::forward<Args>(args)...);
    }
}

template <typename... Args>
void error(spdlog::format_string_t<Args...> fmt, Args&&... args) {
    if (dendrologger) {
        dendrologger->error(fmt, std::forward<Args>(args)...);
    }
}

// "scoped" versions to make the format a bit nicer if we wnt to identify
// internal parts of the program cleanly, this is really only for dendro
// internals
template <typename... Args>
void trace(Scope scope, spdlog::format_string_t<Args...> fmt, Args&&... args) {
    if (dendrologger) {
        // format user message
        auto user_message = fmt::format(fmt, std::forward<Args>(args)...);
        // then log the message
        dendrologger->trace("[{}] {}", scope.name, user_message);
    }
}

template <typename... Args>
void debug(Scope scope, spdlog::format_string_t<Args...> fmt, Args&&... args) {
    if (dendrologger) {
        auto user_message = fmt::format(fmt, std::forward<Args>(args)...);
        dendrologger->debug("[{}] {}", scope.name, user_message);
    }
}

template <typename... Args>
void info(Scope scope, spdlog::format_string_t<Args...> fmt, Args&&... args) {
    if (dendrologger) {
        auto user_message = fmt::format(fmt, std::forward<Args>(args)...);
        dendrologger->info("[{}] {}", scope.name, user_message);
    }
}

template <typename... Args>
void warn(Scope scope, spdlog::format_string_t<Args...> fmt, Args&&... args) {
    if (dendrologger) {
        auto user_message = fmt::format(fmt, std::forward<Args>(args)...);
        dendrologger->warn("[{}] {}", scope.name, user_message);
    }
}

template <typename... Args>
void error(Scope scope, spdlog::format_string_t<Args...> fmt, Args&&... args) {
    if (dendrologger) {
        auto user_message = fmt::format(fmt, std::forward<Args>(args)...);
        dendrologger->error("[{}] {}", scope.name, user_message);
    }
}

}  // namespace logger
}  // namespace dendro
