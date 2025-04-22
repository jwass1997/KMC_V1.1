#include <iostream>
#include <vector>
#include <chrono>

using namespace std;
using namespace std::chrono;

const int ROWS = 5000;
const int COLS = 5000;

// Benchmark using a nested C-style 2D array
void benchmark_array() {
    static int arr[ROWS][COLS];
    // Initialize the array
    for (int i = 0; i < ROWS; ++i)
        for (int j = 0; j < COLS; ++j)
            arr[i][j] = 1;

    volatile long long sum = 0;  // volatile prevents unwanted optimizations
    auto start = high_resolution_clock::now();

    // Sum all elements
    for (int i = 0; i < ROWS; ++i)
        for (int j = 0; j < COLS; ++j)
            sum += arr[i][j];

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start).count();
    cout << "C-style 2D array: Sum = " << sum << ", Time = " << duration << " ms" << endl;
}

// Benchmark using a flattened C-style array (dynamically allocated)
void benchmark_flattened_array() {
    // Allocate one contiguous block for ROWS * COLS integers.
    int* flat_arr = new int[ROWS * COLS];

    // Initialize the flattened array
    for (int i = 0; i < ROWS * COLS; ++i)
        flat_arr[i] = 1;

    volatile long long sum = 0;
    auto start = high_resolution_clock::now();

    // Use manual index arithmetic to simulate 2D access
    for (int i = 0; i < ROWS; ++i)
        for (int j = 0; j < COLS; ++j)
            sum += flat_arr[i * COLS + j];

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start).count();
    cout << "Flattened C-style array: Sum = " << sum << ", Time = " << duration << " ms" << endl;

    // Free the allocated memory
    delete[] flat_arr;
}

// Benchmark using a vector of vectors
void benchmark_nested_vector() {
    // Create a vector of ROWS vectors, each of size COLS, filled with 1
    vector<vector<int>> vec(ROWS, vector<int>(COLS, 1));
    
    volatile long long sum = 0;
    auto start = high_resolution_clock::now();

    // Sum all elements via two levels of indexing
    for (int i = 0; i < ROWS; ++i)
        for (int j = 0; j < COLS; ++j)
            sum += vec[i][j];

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start).count();
    cout << "Nested vector: Sum = " << sum << ", Time = " << duration << " ms" << endl;
}

// Benchmark using a flattened vector with manual index computation
void benchmark_flattened_vector() {
    // Create a single vector with ROWS * COLS elements, filled with 1
    vector<int> vec(ROWS * COLS, 1);
    
    volatile long long sum = 0;
    auto start = high_resolution_clock::now();

    // Access elements using index arithmetic
    for (int i = 0; i < ROWS; ++i)
        for (int j = 0; j < COLS; ++j)
            sum += vec[i * COLS + j];

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start).count();
    cout << "Flattened vector: Sum = " << sum << ", Time = " << duration << " ms" << endl;
}

int main() {
    cout << "Benchmarking a " << ROWS << " x " << COLS << " matrix:" << endl;

    benchmark_array();
    benchmark_flattened_array();
    benchmark_nested_vector();
    benchmark_flattened_vector();

    return 0;
}
