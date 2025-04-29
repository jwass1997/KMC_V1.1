#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>

int main() {

    double min = -1.5;
    double max = 1.5;
    double range = max - min;
    
    int num = 100;

    std::ofstream file;
    file.open("voltages.txt");

    for (int i = 0; i < num; ++i) {
        file << min + i*range/(num-1) << "\n";
    }

    file.close();
}
