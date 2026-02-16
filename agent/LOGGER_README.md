# Logger Class Documentation

A C++20 logging utility for creating component-specific log files with unique identifiers.

## Features

- **Unique Component IDs**: Automatically generates a random 5-digit ID for each logger instance
- **Thread-Safe**: Uses mutex locks for safe multi-threaded logging
- **Multiple Log Levels**: DEBUG, INFO, WARNING, ERROR
- **Automatic Timestamps**: Every log entry includes millisecond-precision timestamps
- **Auto-Flush**: WARNING and ERROR messages are flushed immediately
- **Output Directory**: All logs written to `out/` directory (configurable)
- **Non-Copyable**: Prevents accidental logger duplication

## Basic Usage

### 1. Include the Header

```cpp
#include "logger.hpp"
using namespace logging;
```

### 2. Add Logger as Member Object

```cpp
class MyComponent {
public:
    MyComponent(const std::string& name)
        : logger_(name, LogLevel::INFO)  // Name and optional min level
    {
        logger_.info("MyComponent initialized");
    }

    void do_work() {
        logger_.debug("Starting work...");
        logger_.info("Work completed successfully");
    }

private:
    Logger logger_;  // Logger as member object
};
```

### 3. Log Messages

```cpp
// Basic logging
logger_.debug("Debug message");
logger_.info("Info message");
logger_.warning("Warning message");
logger_.error("Error message");

// Multiple arguments (automatically concatenated)
logger_.info("Processing PC: 0x", std::hex, pc_value);
logger_.warning("Low accuracy: ", accuracy * 100, "%");
logger_.error("Failed with code: ", error_code, " at line ", line_num);
```

## Constructor Parameters

```cpp
Logger(const std::string& component_name,
       LogLevel min_level = LogLevel::INFO,
       const std::string& output_dir = "out")
```

- `component_name`: Name for this component (used in log filename and messages)
- `min_level`: Minimum level to log (default: INFO)
  - `LogLevel::DEBUG`: Log everything
  - `LogLevel::INFO`: Log INFO, WARNING, ERROR
  - `LogLevel::WARNING`: Log WARNING, ERROR only
  - `LogLevel::ERROR`: Log ERROR only
- `output_dir`: Directory for log files (default: "out")

## Log File Format

**Filename**: `<component_name>_<random_id>.log`

**Example**: `TAGE_Predictor_74225.log`

**Content Format**:
```
=== Log started for TAGE_Predictor (ID: 74225) ===
Timestamp: 2026-02-15 15:48:33.399
================================================================================
[2026-02-15 15:48:33.399] [INFO ] [TAGE_Predictor:74225] BranchPredictor initialized
[2026-02-15 15:48:33.400] [DEBUG] [TAGE_Predictor:74225] Component ID: 74225
[2026-02-15 15:48:33.400] [INFO ] [TAGE_Predictor:74225] Branch taken for PC: 0x1234abcd
[2026-02-15 15:48:33.400] [WARN ] [TAGE_Predictor:74225] Low accuracy detected: 75%
```

## Utility Methods

```cpp
// Get component information
std::string name = logger_.get_component_name();
uint32_t id = logger_.get_component_id();
std::string path = logger_.get_log_file_path();

// Change minimum log level at runtime
logger_.set_min_level(LogLevel::DEBUG);
```

## Complete Example

```cpp
#include "logger.hpp"

class TagePredictor {
public:
    TagePredictor()
        : logger_("TAGE", LogLevel::DEBUG)
    {
        logger_.info("Predictor initialized");
        logger_.info("Log file: ", logger_.get_log_file_path());
    }

    void predict(uint64_t pc, bool taken) {
        logger_.debug("Predicting PC: 0x", std::hex, pc);

        // Your prediction logic here...

        if (taken) {
            logger_.info("Branch taken");
        } else {
            logger_.info("Branch not taken");
        }
    }

    void update(bool mispredicted) {
        if (mispredicted) {
            logger_.warning("Misprediction occurred");
            misprediction_count_++;

            if (misprediction_count_ > threshold_) {
                logger_.error("Too many mispredictions: ", misprediction_count_);
            }
        }
    }

private:
    Logger logger_;
    int misprediction_count_ = 0;
    int threshold_ = 1000;
};

int main() {
    TagePredictor predictor;
    predictor.predict(0x1234, true);
    predictor.update(false);
    return 0;
}
```

## Compilation

Requires C++20 and filesystem support:

```bash
g++ -std=c++20 -o myprogram myprogram.cpp
```

## Thread Safety

The logger is thread-safe. Multiple threads can call logging methods on the same logger instance without data races.

## Performance Notes

- DEBUG and INFO messages are buffered (flushed periodically)
- WARNING and ERROR messages are flushed immediately
- Consider using DEBUG level only during development (set to INFO for production)

## Best Practices

1. **One Logger Per Component**: Each major component should have its own logger
2. **Descriptive Names**: Use clear component names (e.g., "L1_Cache", "TAGE_Predictor")
3. **Appropriate Levels**:
   - DEBUG: Detailed flow, variable values
   - INFO: Normal operations, state changes
   - WARNING: Unexpected but handled conditions
   - ERROR: Failures, critical issues
4. **Structured Messages**: Include context (PC values, IDs, etc.)
5. **Set Min Level**: Use `LogLevel::INFO` for production, `DEBUG` for debugging

## Integration with Existing Code

For your TageTable class:

```cpp
class TageTable {
public:
    TageTable()
        : logger_("TageTable", LogLevel::INFO)
    {
        logger_.info("TageTable initialized with ", num_entries_, " entries");
    }

    void lookup(uint64_t pc) {
        logger_.debug("Lookup for PC: 0x", std::hex, pc);
        // ... lookup logic
    }

private:
    Logger logger_;
    size_t num_entries_;
};
```
