#include <vector>
#include <iostream>
#include <tuple>
#include <chrono>
#include <string>

#include "Simulator.h"
#include "utils.h"
#include "SystemGraph.h"

Simulator::Simulator() 
    : system(nullptr) 
{
    system = new SystemGraph();
}

Simulator::Simulator(const std::string& path) 
    : system(nullptr) 
{
    if(path.empty()) {
        throw std::invalid_argument("Empty path");
    }

    system = new SystemGraph(path);
}

Simulator::Simulator(const std::string& path, int burnInHops, int numOfHops) 
    : system(nullptr) 
    , burnInHops(burnInHops)
    , numOfHops(numOfHops)
{
    if(path.empty()) {
        throw std::invalid_argument("Empty path");
    }
    
    system = new SystemGraph(path);
}

Simulator::~Simulator() {
    delete system;
    system = nullptr;
}

void Simulator::singleMCStep(bool trackingData) {

    int numOfStates = system->getNumOfStates();

    for (int indexOfState = 0; indexOfState < numOfStates; ++indexOfState) {
        system->updateStateEnergy(indexOfState);
    }

    system->updateTransitionRates();

    auto [selectedI, selectedJ, totalRate] = system->sampleTransitionEvent();

    system->updateStateOccupation(selectedI, selectedJ);

    system->increaseSystemTime(totalRate);

    if(trackingData) {
        system->eventCounts[selectedI*numOfStates + selectedJ] += 1;
    }
}

void Simulator::simulateNumberOfSteps(unsigned int numOfSteps, bool trackingData) {

    /**
     * 
     * For equilibrium simulation set trackingData = false
     * 
     */

    for(int i = 0; i < numOfSteps; ++i) {
        singleMCStep(trackingData);
    }
}