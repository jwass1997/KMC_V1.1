#include <iostream>
#include <vector>
#include <filesystem>
#include <string>
#include <chrono>
#include <omp.h>

#include "SystemGraph.h"
#include "Simulator.h"
#include "utils.h"

int main(int argc, char* argv[]) {
    std::string saveFolder = "currentData";
    std::string defaultDeviceConfigs = "default_configs";
    createDirectoryFromStringPath("", saveFolder);
    int controlElectrodeIndex = 1;
    int scanElectrodeIndex = 0;

    std::string deviceID1 = "1";
    std::vector<int> inputElectrodes = {1,2,3,4,5,6,7};
    std::vector<int> outpuElectrodes = {0};
    int maxThreads = omp_get_max_threads();
    //std::cout << maxThreads << "\n";
    //std::cout<<__cplusplus<<"\n";
    //createBatchOfSingleSystem(20, inputElectrodes, outpuElectrodes, -1.5, 1.5, 1e4, 1e6, 1, defaultDeviceConfigs, saveFolder, 1);
    //recordDevice(deviceID1, 1e4, 1e5, defaultDeviceConfigs, saveFolder);

    int numOfPoints = 100;
    double maxVoltage = -1.;
    double minVoltage = 1.;
    double range = maxVoltage - minVoltage;
    double voltageStep = range / (numOfPoints-1);
    std::vector<double> currents(numOfPoints, 0.0);
    #pragma omp parallel for
    for (int voltagePoint = 0; voltagePoint < numOfPoints; ++voltagePoint) {
        std::cout<<voltagePoint<<"\n";
        std::vector<double> voltageSetting = {0.,0.,0.,0.,0.,0.,0.,0.};
        voltageSetting[1] = minVoltage + voltageStep*voltagePoint;
        double outputCurrent = 0.0;
        outputCurrent = IVPoint(voltageSetting, 5, 0, 1e4, 1e7, 100, voltagePoint, defaultDeviceConfigs, saveFolder);
        currents[voltagePoint] = outputCurrent;
    }
    std::string fileName = "currents.npz";
    cnpy::npz_save(fileName, "currents", currents.data(), {static_cast<size_t>(numOfPoints)}, "w");
} 