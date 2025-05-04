#include <iostream>
#include <vector>
#include <filesystem>
#include <string>
#include <chrono>
#include <omp.h>
#include <fstream>

#include "SystemGraph.h"
#include "Simulator.h"
#include "FEMmethods.h"
#include "utils.h"

int main(int argc, char* argv[]) {

    std::string defaultConfig = "default_configs";

    /* int numPoints = 100;
    double minVoltage = -1.5;
    double maxVoltage = 1.5;
    double range = maxVoltage - minVoltage;
    double vStep = range / static_cast<double>(numPoints-1);    

    int numOfDevices = 5;
    int scanElectrodeIndex = 0;
    int controlElectrode = 1;
    int equilibriumSteps = 1e4;
    int simulationSteps = 1e5;
    int numOfIntervals = 1;
    
    std::ofstream file;
    file.open("iv_curve.csv");
    double current = 0.0;
    std::vector<double> voltageSetting = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (int _point = 0; _point < numPoints; ++_point) {

        voltageSetting[controlElectrode] = minVoltage + -_point*vStep;
        current = IVPoint(
            voltageSetting,
            numOfDevices,
            scanElectrodeIndex,
            equilibriumSteps,
            simulationSteps,
            numOfIntervals,
            defaultConfig  
        );

        std::cout << "_point==" << _point << "\n"
                  << "current==" << current << "\n";

        file << voltageSetting[controlElectrode] << current << "\n";
    }

    file.close(); */

    argParser(argc, argv);
} 