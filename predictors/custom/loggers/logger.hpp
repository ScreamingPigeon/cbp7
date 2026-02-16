#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <fstream>
#include <string>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <random>
#include <mutex>
#include <filesystem>
#include <ctime>

namespace logging {

enum class LogLevel {
    DEBUG,
    INFO,
    WARNING,
    ERROR
};

class Logger {
public:
    // Constructor: takes a component name
    explicit Logger(const std::string& component_name,
                   LogLevel min_level = LogLevel::INFO,
                   const std::string& output_dir = "out")
        : component_name_(component_name)
        , component_id_(generate_random_id())
        , min_level_(min_level)
        , output_dir_(output_dir)
    {
        initialize_log_file();
    }

    // Destructor: ensure all logs are flushed
    ~Logger() {
        if (log_file_.is_open()) {
            log_file_.flush();
            log_file_.close();
        }
    }

    // Delete copy and move operations (prevent copying and moving due to mutex)
    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;
    Logger(Logger&&) = delete;
    Logger& operator=(Logger&&) = delete;

    // Logging methods
    template<typename... Args>
    void debug(Args&&... args) {
        log(LogLevel::DEBUG, std::forward<Args>(args)...);
    }

    template<typename... Args>
    void info(Args&&... args) {
        log(LogLevel::INFO, std::forward<Args>(args)...);
    }

    template<typename... Args>
    void warning(Args&&... args) {
        log(LogLevel::WARNING, std::forward<Args>(args)...);
    }

    template<typename... Args>
    void error(Args&&... args) {
        log(LogLevel::ERROR, std::forward<Args>(args)...);
    }

    // Get component info
    std::string get_component_name() const { return component_name_; }
    uint32_t get_component_id() const { return component_id_; }
    std::string get_log_file_path() const { return log_file_path_; }

    // Set minimum log level
    void set_min_level(LogLevel level) { min_level_ = level; }

private:
    std::string component_name_;
    uint32_t component_id_;
    LogLevel min_level_;
    std::string output_dir_;
    std::string log_file_path_;
    std::ofstream log_file_;
    std::mutex mutex_;  // For thread-safety

    // Generate a random component ID
    uint32_t generate_random_id() {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<uint32_t> dis(10000, 99999);
        return dis(gen);
    }

    // Initialize log file
    void initialize_log_file() {
        // Create output directory if it doesn't exist
        std::filesystem::create_directories(output_dir_);

        // Generate log file name: <component_name>_<component_id>.log
        std::ostringstream filename;
        filename << output_dir_ << "/" << component_name_ << "_" << component_id_ << ".log";
        log_file_path_ = filename.str();

        // Open log file
        log_file_.open(log_file_path_, std::ios::out | std::ios::app);

        if (!log_file_.is_open()) {
            throw std::runtime_error("Failed to open log file: " + log_file_path_);
        }

        // Write header
        log_file_ << "=== Log started for " << component_name_
                  << " (ID: " << component_id_ << ") ===" << std::endl;
        log_file_ << "Timestamp: " << get_timestamp() << std::endl;
        log_file_ << std::string(80, '=') << std::endl;
        log_file_.flush();
    }

    // Get current timestamp as string
    std::string get_timestamp() {
        auto now = std::chrono::system_clock::now();
        auto now_time_t = std::chrono::system_clock::to_time_t(now);
        auto now_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            now.time_since_epoch()) % 1000;

        std::ostringstream oss;
        oss << std::put_time(std::localtime(&now_time_t), "%Y-%m-%d %H:%M:%S");
        oss << '.' << std::setfill('0') << std::setw(3) << now_ms.count();
        return oss.str();
    }

    // Get log level string
    const char* level_to_string(LogLevel level) {
        switch (level) {
            case LogLevel::DEBUG:   return "DEBUG";
            case LogLevel::INFO:    return "INFO ";
            case LogLevel::WARNING: return "WARN ";
            case LogLevel::ERROR:   return "ERROR";
            default:                return "UNKNOWN";
        }
    }

    // Main logging function
    template<typename... Args>
    void log(LogLevel level, Args&&... args) {
        // Check if this level should be logged
        if (level < min_level_) {
            return;
        }

        // Thread-safe logging
        std::lock_guard<std::mutex> lock(mutex_);

        if (!log_file_.is_open()) {
            return;
        }

        // Format: [TIMESTAMP] [LEVEL] [COMPONENT_NAME:ID] Message
        log_file_ << "[" << get_timestamp() << "] "
                  << "[" << level_to_string(level) << "] "
                  << "[" << component_name_ << ":" << component_id_ << "] ";

        // Write all arguments
        ((log_file_ << std::forward<Args>(args)), ...);
        log_file_ << std::endl;

        // Flush immediately for important messages
        if (level >= LogLevel::WARNING) {
            log_file_.flush();
        }
    }
};

} // namespace logging

#endif // LOGGER_HPP
