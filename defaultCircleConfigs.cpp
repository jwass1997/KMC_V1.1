#include <vector>
#include <iostream>
#include <limits>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <iostream>
#include <random>

std::random_device rd;
std::mt19937 rng(rd());

std::uniform_real_distribution<double> uniformDistribution(0.0, 1.0);

int main() {
    std::cout<<__cplusplus<<"\n";
    int nAcceptors = 200;
    int nDonors = 3;
    int nElectrodes = 8;
    double radius = 100.0;
    double R = std::sqrt(M_PI*radius*radius / static_cast<double>(nAcceptors));
    
    radius = radius / R;

    std::vector<double> defaultAngularPositions = {
        0.0,
        45.0,
        90.0,
        135.0,
        180.0,
        225.0,
        270.0,
        315.0
    };

    std::vector<double> defaultVoltages = {
        1.5,
        1.2,
        1.1,
        0.9,
        1.2,
        -0.8,
        -1.2,
        -1.0
    };

    std::ofstream file;
    file.open("default_configs/default_circle_acceptors.txt");

    for(int i = 0; i < nAcceptors; ++i) {
        double angle = 2.0*M_PI*uniformDistribution(rng);
        double r = radius*std::sqrt(uniformDistribution(rng));

        double xCoordinate = r*std::cos(angle);
        double yCoordinate = r*std::sin(angle);

        file << xCoordinate << " " << yCoordinate << "\n";
    }
    file.close();

    file.open("default_configs/default_circle_donors.txt");
    for(int i = 0; i < nDonors; ++i) {
        double angle = 2.0*M_PI*uniformDistribution(rng);
        double r = radius*std::sqrt(uniformDistribution(rng));

        double xCoordinate = r*std::cos(angle);
        double yCoordinate = r*std::sin(angle);

        file << xCoordinate << " " << yCoordinate << "\n";
    }
    file.close();

    file.open("default_configs/default_circle_electrodes.txt");
    for(int i = 0; i < nElectrodes; ++i) {
        double defaultVoltage = 1.5*uniformDistribution(rng) - 1.5;
        file << defaultAngularPositions[i] << " " << defaultVoltage << "\n";
    }
    file.close();
}