#include "logger.hpp"
#include <iostream>
#include <thread>
#include <chrono>

using namespace logging;

// Example class 1: Branch predictor component
class BranchPredictor {
public:
    BranchPredictor(const std::string& name)
        : logger_(name, LogLevel::DEBUG)  // Initialize logger with component name
    {
        logger_.info("BranchPredictor initialized");
        logger_.debug("Component ID: ", logger_.get_component_id());
    }

    void predict(uint64_t pc, bool taken) {
        logger_.debug("Predicting for PC: 0x", std::hex, pc, std::dec);

        // Simulate prediction
        if (taken) {
            logger_.info("Branch taken for PC: 0x", std::hex, pc);
        } else {
            logger_.info("Branch not taken for PC: 0x", std::hex, pc);
        }
    }

    void report_accuracy(double accuracy) {
        if (accuracy < 0.8) {
            logger_.warning("Low accuracy detected: ", accuracy * 100, "%");
        } else {
            logger_.info("Accuracy: ", accuracy * 100, "%");
        }
    }

    void handle_error(const std::string& error_msg) {
        logger_.error("Error occurred: ", error_msg);
    }

private:
    Logger logger_;  // Logger as member object
};

// Example class 2: Cache simulator
class CacheSimulator {
public:
    CacheSimulator(const std::string& cache_type)
        : logger_(cache_type + "_Cache", LogLevel::INFO)
    {
        logger_.info("Cache simulator created for ", cache_type);
        logger_.info("Log file: ", logger_.get_log_file_path());
    }

    void access(uint64_t address, bool hit) {
        if (hit) {
            logger_.debug("Cache hit at address: 0x", std::hex, address);
        } else {
            logger_.warning("Cache miss at address: 0x", std::hex, address);
        }
    }

    void flush() {
        logger_.info("Cache flushed");
    }

private:
    Logger logger_;
};

// Example class 3: Performance monitor
class PerformanceMonitor {
public:
    PerformanceMonitor()
        : logger_("PerfMonitor", LogLevel::DEBUG)
    {
        logger_.info("Performance monitoring started");
        start_time_ = std::chrono::steady_clock::now();
    }

    void log_metric(const std::string& metric_name, double value) {
        logger_.info("Metric [", metric_name, "]: ", value);
    }

    void checkpoint(const std::string& checkpoint_name) {
        auto now = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
            now - start_time_).count();
        logger_.info("Checkpoint [", checkpoint_name, "] at ", elapsed, "ms");
    }

private:
    Logger logger_;
    std::chrono::steady_clock::time_point start_time_;
};

int main() {
    std::cout << "Logger Example - Creating multiple components with loggers\n";
    std::cout << "Check the 'out/' directory for log files\n\n";

    // Create instances of different classes, each with their own logger
    BranchPredictor bp1("TAGE_Predictor");
    BranchPredictor bp2("Bimodal_Predictor");
    CacheSimulator l1_cache("L1");
    CacheSimulator l2_cache("L2");
    PerformanceMonitor perf_mon;

    // Each component logs independently
    std::cout << "\n=== Branch Predictor 1 ===\n";
    bp1.predict(0x1234ABCD, true);
    bp1.predict(0x5678EF01, false);
    bp1.report_accuracy(0.85);

    std::cout << "\n=== Branch Predictor 2 ===\n";
    bp2.predict(0xDEADBEEF, true);
    bp2.report_accuracy(0.75);
    bp2.handle_error("Prediction buffer overflow");

    std::cout << "\n=== Cache Simulators ===\n";
    l1_cache.access(0x1000, true);
    l1_cache.access(0x2000, false);
    l2_cache.access(0x3000, false);
    l2_cache.flush();

    std::cout << "\n=== Performance Monitor ===\n";
    perf_mon.log_metric("IPC", 2.5);
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    perf_mon.checkpoint("After initialization");
    perf_mon.log_metric("Cache miss rate", 0.15);

    std::cout << "\n=== Done ===\n";
    std::cout << "All logs written to 'out/' directory\n";
    std::cout << "Each component has a unique log file with its random ID\n";

    return 0;
}
