#include <iostream>
#include <vector>
#include <filesystem>
#include <string>
#include <chrono>

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
    std::cout<<__cplusplus<<"\n";
    //createBatchOfSingleSystem(20, inputElectrodes, outpuElectrodes, -1.5, 1.5, 1e4, 1e4, 1, defaultDeviceConfigs, saveFolder, 1);
    //auto start = std::chrono::high_resolution_clock::now();
    //recordDevice(deviceID1, 1e4, 1e5, defaultDeviceConfigs, saveFolder);
    //auto end = std::chrono::high_resolution_clock::now();
    //std::cout << "Simulation time: " << std::chrono::duration_cast<std::chrono::seconds>(end-start).count() << "\n";
    voltageCurrentCharacteristic(-1.5, 1.5, 20, controlElectrodeIndex, scanElectrodeIndex, 1e4, 1e5, 1, defaultDeviceConfigs, saveFolder);
}