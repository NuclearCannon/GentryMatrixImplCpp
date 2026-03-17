#pragma once
#include <chrono>

class HighResolutionTimer {
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
    
public:
    void start() {
        start_time = std::chrono::high_resolution_clock::now();
    }
    
    double stop() {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
            end_time - start_time
        );
        return duration.count() / 1000000.0; // 返回毫秒
    }
};

