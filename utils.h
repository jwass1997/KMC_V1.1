#pragma once

#include <iostream>
#include <random>
#include <string>
#include <filesystem>

#include "cnpy.h"

class SystemGraph;
class Simulator;

struct Electrode {
    double angularPosition;
    double voltage;
};

extern std::mt19937 rng;

inline constexpr double kb = 1.380649e-23;

inline constexpr double e = 1.602176634e-19;

inline constexpr double PI = 3.1415926535897;

inline constexpr double eps0 = 8.854187817620389e-12;

inline constexpr double epsr = 10.0;

inline constexpr double kbT = 1.0630997e-21;

double sampleFromUniformDistribution(double min, double max);

double sampleFromNormalDistribution(double mean, double standardDeviation);

double calculateDistance(double coordinateX1, double coordinateX2, double coordinateY1, double coordinateY2);

void createDirectoryFromStringPath(const std::string& path, const std::string& directoryName);

void recordDevice(const std::string& ID, int equilibriumSteps, int numOfSteps, const std::string& defaultConfigs, const std::string& saveFolderPath);

void voltageCurrentCharacteristic(
    double minVoltage, 
    double maxVoltage, 
    int numOfPoints, 
    int controlElectrodeIndex, 
    int scanElectrodeIndex, 
    int equilibriumSteps, 
    int scanSteps, 
    int ID, 
    const std::string& defaultConfigs, 
    const std::string& saveFolderPath
);

double calculateCurrentAverage(
    Simulator& simulator, 
    int controlElectrodeIndex, 
    int scanElectrodeIndex, 
    double voltage, 
    int equilibriumSteps, 
    int simulationSteps, 
    int numberOfIntervals
);

void createBatchOfSingleSystem(
    int batchSize, 
    std::vector<int> inputElectrodes,
    std::vector<int> outputElectrodes,
    double minVoltage,
    double maxVoltage,
    int equilibriumSteps, 
    int simulationSteps, 
    int numOfIntervals,
    const std::string& defaultConfigs, 
    const std::string& saveFolder, 
    int batchID
);