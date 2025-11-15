#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <deque>
#include <list>
#include <chrono>
#include <random>
#include <numeric>
#include <iomanip>
#include <stdexcept>

enum class ContainerType { Vector, Deque, List };
enum class Operation {
    PushBack, PushFront, RandomInsert,
    RandomErase, IterateSum, RandomAccess
};
enum class TimeUnit { Milliseconds, Seconds };

union TimeValue {
    double ms;
    double sec;
};

struct TimeResult {
    Operation op;
    ContainerType ct;
    int N;
    TimeUnit unit;
    TimeValue val;
};

std::string op_to_string(Operation op) {
    switch (op) {
        case Operation::PushBack:     return "PushBack";
        case Operation::PushFront:    return "PushFront";
        case Operation::RandomInsert: return "RandomInsert";
        case Operation::RandomErase:  return "RandomErase";
        case Operation::IterateSum:   return "IterateSum";
        case Operation::RandomAccess: return "RandomAccess";
        default:                      return "Unknown";
    }
}

std::string ct_to_string(ContainerType ct) {
    switch (ct) {
        case ContainerType::Vector: return "Vector";
        case ContainerType::Deque:  return "Deque";
        case ContainerType::List:   return "List";
        default:                    return "Unknown";
    }
}

std::vector<int> generate_payload(int N) {
    std::vector<int> payload(N);
    std::mt19937 gen(12345);
    std::uniform_int_distribution<int> dist(0, 1000000);
    for (int& x : payload) {
        x = dist(gen);
    }
    return payload;
}

std::vector<int> generate_indices(int N, int max_index) {
    if (max_index == 0) return std::vector<int>(N, 0);
    std::vector<int> indices(N);
    std::mt19937 gen(54321);
    std::uniform_int_distribution<int> dist(0, max_index - 1);
    for (int& x : indices) {
        x = dist(gen);
    }
    return indices;
}

template<class Seq>
int access_element(const Seq& c, int index) {
    return c[index];
}

int access_element(const std::list<int>& c, int index) {
    return *std::next(c.begin(), index);
}

template<class Seq>
auto get_iterator(Seq& c, int index) -> decltype(c.begin()) {
    return c.begin() + index;
}

auto get_iterator(std::list<int>& c, int index) -> decltype(c.begin()) {
    return std::next(c.begin(), index);
}

template<class Seq>
double measure_push_back(const std::vector<int>& payload, int repeats) {
    double total_time = 0.0;
    for (int i = 0; i < repeats; ++i) {
        Seq c;
        auto start = std::chrono::high_resolution_clock::now();
        for (int x : payload) {
            c.push_back(x);
        }
        auto end = std::chrono::high_resolution_clock::now();
        total_time += std::chrono::duration<double, std::milli>(end - start).count();
    }
    return total_time / repeats;
}

template<class Seq>
double measure_push_front(const std::vector<int>& payload, int repeats) {
    double total_time = 0.0;
    for (int i = 0; i < repeats; ++i) {
        Seq c;
        auto start = std::chrono::high_resolution_clock::now();
        for (int x : payload) {
            c.push_front(x);
        }
        auto end = std::chrono::high_resolution_clock::now();
        total_time += std::chrono::duration<double, std::milli>(end - start).count();
    }
    return total_time / repeats;
}

double measure_push_front_vector(const std::vector<int>& payload, int repeats) {
    double total_time = 0.0;
    for (int i = 0; i < repeats; ++i) {
        std::vector<int> c;
        auto start = std::chrono::high_resolution_clock::now();
        for (int x : payload) {
            c.insert(c.begin(), x);
        }
        auto end = std::chrono::high_resolution_clock::now();
        total_time += std::chrono::duration<double, std::milli>(end - start).count();
    }
    return total_time / repeats;
}

template<class Seq>
double measure_random_insert(const std::vector<int>& payload, int repeats) {
    int N = payload.size();
    double total_time = 0.0;
    std::mt19937 gen(12345);

    for (int i = 0; i < repeats; ++i) {
        Seq c;
        std::uniform_int_distribution<int> dist(0, 0);
        using dist_params = typename std::uniform_int_distribution<int>::param_type;

        auto start = std::chrono::high_resolution_clock::now();
        for (int j = 0; j < N; ++j) {
            int pos = (c.empty()) ? 0 : dist(gen, dist_params(0, c.size()));
            auto it = get_iterator(c, pos);
            c.insert(it, payload[j]);
        }
        auto end = std::chrono::high_resolution_clock::now();
        total_time += std::chrono::duration<double, std::milli>(end - start).count();
    }
    return total_time / repeats;
}


template<class Seq>
double measure_random_erase(const std::vector<int>& payload, int repeats) {
    int N = payload.size();
    double total_time = 0.0;
    std::mt19937 gen(12345);

    for (int i = 0; i < repeats; ++i) {
        Seq c(payload.begin(), payload.end());
        std::uniform_int_distribution<int> dist(0, 0);
        using dist_params = typename std::uniform_int_distribution<int>::param_type;

        auto start = std::chrono::high_resolution_clock::now();
        for (int j = 0; j < N; ++j) {
            int pos = dist(gen, dist_params(0, c.size() - 1));
            auto it = get_iterator(c, pos);
            c.erase(it);
        }
        auto end = std::chrono::high_resolution_clock::now();
        total_time += std::chrono::duration<double, std::milli>(end - start).count();
    }
    return total_time / repeats;
}


template<class Seq>
double measure_iterate_sum(const std::vector<int>& payload, int repeats) {
    double total_time = 0.0;
    Seq c(payload.begin(), payload.end());

    for (int i = 0; i < repeats; ++i) {
        auto start = std::chrono::high_resolution_clock::now();

        volatile long long sum = 0;
        for (int x : c) {
            sum += x;
        }

        auto end = std::chrono::high_resolution_clock::now();
        total_time += std::chrono::duration<double, std::milli>(end - start).count();
    }
    return total_time / repeats;
}


template<class Seq>
double measure_random_access(const std::vector<int>& payload, int repeats) {
    int N = payload.size();
    double total_time = 0.0;
    Seq c(payload.begin(), payload.end());
    std::vector<int> indices = generate_indices(N, N);

    for (int i = 0; i < repeats; ++i) {
        auto start = std::chrono::high_resolution_clock::now();

        volatile int val = 0;
        for (int idx : indices) {
            val = access_element(c, idx);
        }

        auto end = std::chrono::high_resolution_clock::now();
        total_time += std::chrono::duration<double, std::milli>(end - start).count();
    }
    return total_time / repeats;
}

int main() {
    const std::vector<int> sizes = { 8192, 32768, 262144, 500000 };
    const int base_repeats = 10;

    std::vector<TimeResult> all_results;

    std::ofstream result_file("result.txt");
    if (!result_file.is_open()) {
        std::cerr << "Error\n"; 
        return 1;
    }

    std::cout << std::setw(14) << std::left << "Operation";
    std::cout << std::setw(10) << std::left << "Container";
    std::cout << std::setw(12) << std::right << "N";
    std::cout << std::setw(16) << std::right << "Time (ms)";
    std::cout << '\n'; 
    std::cout << std::string(52, '-') << '\n'; 
    result_file << std::setw(14) << std::left << "Operation";
    result_file << std::setw(10) << std::left << "Container";
    result_file << std::setw(12) << std::right << "N";
    result_file << std::setw(16) << std::right << "Time (ms)";
    result_file << '\n'; 
    result_file << std::string(52, '-') << '\n'; 
    result_file << std::fixed << std::setprecision(6);
    std::cout << std::fixed << std::setprecision(6);

    for (int N : sizes) {
        std::cout << "N=" << N << '\n'; 
        int slow_repeats = (N >= 262144) ? 1 : base_repeats;
        std::vector<int> payload = generate_payload(N);

        all_results.push_back({Operation::PushBack, ContainerType::Vector, N, TimeUnit::Milliseconds, {measure_push_back<std::vector<int>>(payload, base_repeats)}});
        all_results.push_back({Operation::PushBack, ContainerType::Deque,  N, TimeUnit::Milliseconds, {measure_push_back<std::deque<int>>(payload, base_repeats)}});
        all_results.push_back({Operation::PushBack, ContainerType::List,   N, TimeUnit::Milliseconds, {measure_push_back<std::list<int>>(payload, base_repeats)}});

        all_results.push_back({Operation::PushFront, ContainerType::Vector, N, TimeUnit::Milliseconds, {measure_push_front_vector(payload, slow_repeats)}});
        all_results.push_back({Operation::PushFront, ContainerType::Deque,  N, TimeUnit::Milliseconds, {measure_push_front<std::deque<int>>(payload, base_repeats)}});
        all_results.push_back({Operation::PushFront, ContainerType::List,   N, TimeUnit::Milliseconds, {measure_push_front<std::list<int>>(payload, base_repeats)}});

        all_results.push_back({Operation::RandomInsert, ContainerType::Vector, N, TimeUnit::Milliseconds, {measure_random_insert<std::vector<int>>(payload, slow_repeats)}});
        all_results.push_back({Operation::RandomInsert, ContainerType::Deque,  N, TimeUnit::Milliseconds, {measure_random_insert<std::deque<int>>(payload, slow_repeats)}});
        all_results.push_back({Operation::RandomInsert, ContainerType::List,   N, TimeUnit::Milliseconds, {measure_random_insert<std::list<int>>(payload, slow_repeats)}});

        all_results.push_back({Operation::RandomErase, ContainerType::Vector, N, TimeUnit::Milliseconds, {measure_random_erase<std::vector<int>>(payload, slow_repeats)}});
        all_results.push_back({Operation::RandomErase, ContainerType::Deque,  N, TimeUnit::Milliseconds, {measure_random_erase<std::deque<int>>(payload, slow_repeats)}});
        all_results.push_back({Operation::RandomErase, ContainerType::List,   N, TimeUnit::Milliseconds, {measure_random_erase<std::list<int>>(payload, slow_repeats)}});

        all_results.push_back({Operation::IterateSum, ContainerType::Vector, N, TimeUnit::Milliseconds, {measure_iterate_sum<std::vector<int>>(payload, base_repeats)}});
        all_results.push_back({Operation::IterateSum, ContainerType::Deque,  N, TimeUnit::Milliseconds, {measure_iterate_sum<std::deque<int>>(payload, base_repeats)}});
        all_results.push_back({Operation::IterateSum, ContainerType::List,   N, TimeUnit::Milliseconds, {measure_iterate_sum<std::list<int>>(payload, base_repeats)}});

        all_results.push_back({Operation::RandomAccess, ContainerType::Vector, N, TimeUnit::Milliseconds, {measure_random_access<std::vector<int>>(payload, base_repeats)}});
        all_results.push_back({Operation::RandomAccess, ContainerType::Deque,  N, TimeUnit::Milliseconds, {measure_random_access<std::deque<int>>(payload, base_repeats)}});
        all_results.push_back({Operation::RandomAccess, ContainerType::List,   N, TimeUnit::Milliseconds, {measure_random_access<std::list<int>>(payload, slow_repeats)}});

        for (size_t i = all_results.size() - 18; i < all_results.size(); ++i) {
            TimeResult& r = all_results[i];
            std::cout << std::setw(14) << std::left << op_to_string(r.op);
            std::cout << std::setw(10) << std::left << ct_to_string(r.ct);
            std::cout << std::setw(12) << std::right << r.N;
            std::cout << std::setw(16) << std::right << r.val.ms;
            std::cout << '\n'; 
            result_file << std::setw(14) << std::left << op_to_string(r.op);
            result_file << std::setw(10) << std::left << ct_to_string(r.ct);
            result_file << std::setw(12) << std::right << r.N;
            result_file << std::setw(16) << std::right << r.val.ms;
            result_file << '\n'; 
        }
    }
    result_file.close();
    return 0;
}