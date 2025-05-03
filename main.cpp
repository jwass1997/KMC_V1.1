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
    /* //argParser(argc, argv);
    

    int nAcceptors = 200;
    int numOfElectrodes = 8;
    int numOfStates = 200 + numOfElectrodes;
    int scanElectrodeIndex = 0;
    int simulationSteps = 1e5;
    int numOfIntervals = 100;
    Simulator sim(defaultConfig);

    sim.simulateNumberOfSteps(1e4, false);

    double averageCurrent = 0.0;
    double totalTime = 0.0;
    int totalNet = 0;
    int intervalSteps = simulationSteps / numOfIntervals;
    int intervalCounter = 0;
    while(intervalCounter < numOfIntervals) {
        double startClock = sim.system->getSystemTime();
        sim.simulateNumberOfSteps(intervalSteps, true);
        double endClock = sim.system->getSystemTime();
        double elapsedTime = endClock - startClock;
        int inCounts = 0;
        int outCounts = 0;
        for (int i = 0; i < numOfStates; ++i) {
            outCounts += sim.system->getNumberOfEvents(nAcceptors+scanElectrodeIndex, i);
            inCounts += sim.system->getNumberOfEvents(i, nAcceptors+scanElectrodeIndex); 
        }
        totalTime += elapsedTime;
        totalNet += inCounts-outCounts;

        sim.system->resetEventCounts();
        intervalCounter++;
    }
    averageCurrent = static_cast<double>(totalNet) / totalTime;

    std::cout<<averageCurrent<<"\n"; */
    Simulator sim(defaultConfig);

    int nA = sim.system->nAcceptors;
    int nE = sim.system->nElectrodes;
    int numOfStates = 200 + nE;
    int numPoints = 20;
    double minVoltage = -1.5;
    double maxVoltage = 1.5;
    double range = maxVoltage - minVoltage;
    double vStep = range / static_cast<double>(numPoints-1);

    double current = 0.0;
    std::ofstream file;
    file.open("current_save.csv");
    std::vector<double> voltageSetting = {
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0
    };
    /* int outputElectrodeIdx = 0;
    for (int i = 0; i < numPoints; ++i) {
        voltageSetting[1] = minVoltage + i*vStep;
        //std::cout << voltageSetting[1] << "\n";
        current = 0.0;
        sim.simulateNumberOfSteps(1e4, false);
        sim.system->updateVoltages(voltageSetting);
        double startTime = sim.system->getSystemTime();
        sim.simulateNumberOfSteps(1e5, true);
        int inCounts = 0;
        int outCounts = 0;
        for (int i = 0; i < numOfStates; ++i) {
            inCounts += sim.system->getNumberOfEvents(i, nA+outputElectrodeIdx);
            outCounts += sim.system->getNumberOfEvents(nA+outputElectrodeIdx, i);
        }
        current = static_cast<double>(inCounts - outCounts);
        double endTime = sim.system->getSystemTime();
        file << current / (endTime - startTime) << "\n";
    }
    file.close(); */
    //recordDevice("1", 1e4, 1e6, defaultConfig, "currentData");
    /* for (int i = 0; i < sim.system->numOfStates; ++i) {
        std::cout << sim.system->stateEnergies[i] << "\n";
    } */
    //sim.simulateNumberOfSteps(1e3, true);
    //std::vector<int> netCounts(numOfElectrodes, 0);
    /* for (int i = 0; i < numOfElectrodes; ++i) {
        int inCounts = 0;
        int outCounts = 0;
        for (int j = 0; j < numOfStates; ++j) {
            inCounts += sim.system->getNumberOfEvents(j, nAcceptors+i);
            outCounts += sim.system->getNumberOfEvents(nAcceptors+i, j);
        }
        netCounts[i] = inCounts - outCounts;
        std::cout << netCounts[i] << "\n";
    } */

    int nS = sim.system->numOfStates;

    /* for (int i = 0; i < nS; ++i) {
        std::cout << sim.system->stateEnergies[i] << "\n";
    } */

    /* auto fem = sim.system->finiteElementSolver;
    auto mesh = fem->mesh;
    int numV = mesh->GetNV();
    std::ofstream out("my_potential.csv");
    out << "x,y,potential\n";

    for (int i = 0; i < numV; ++i) {
        const auto &v = mesh->GetVertex(i);
        double x = v[0];
        double y = v[1];
        double phi = fem->calculatePotential(x, y);
        out << x << "," << y << "," << phi << "\n";
    }

    out.close(); */

    /* for (int i = 0; i < nAcceptors; ++i) {
        std::cout << sim.system->acceptorCoordinates[i*2] << " " << sim.system->acceptorCoordinates[i*2+1] << "\n";
    } */
} 