#include <iostream>
#include <vector>
#include <filesystem>
#include <string>
#include <chrono>
#include <omp.h>

#include "SystemGraph.h"
#include "Simulator.h"
#include "utils.h"

int main() {
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
    recordDevice(deviceID1, 1e4, 1e6, defaultDeviceConfigs, saveFolder);

    //IVCurve(-1.5, 1.5, 1, 100, 0, 1e4, 1e6, 10, 1, defaultDeviceConfigs, saveFolder);
    /* int numOfPoints = 2;
    std::vector<double> currents(numOfPoints, 0.0);
    double maxVoltage = -1.5;
    double minVoltage = 1.5;
    double range = maxVoltage - minVoltage;
    double voltageStep = range / (numOfPoints-1);
    std::vector<double> voltageSetting = sampleVoltageSetting(8, minVoltage, maxVoltage);
    //#pragma omp parallel for
    for (int i = 0; i < numOfPoints; ++i) {
        voltageSetting[1] = minVoltage + voltageStep*i;
        double current = currentFromVoltageCombination(voltageSetting, 0, 1e4, 1e6, 10, defaultDeviceConfigs);
        currents[i] = current;
    }
    std::string fileName = "currents.npz";
    cnpy::npz_save(fileName, "currents", currents.data(), {static_cast<size_t>(numOfPoints)}, "w"); */
}