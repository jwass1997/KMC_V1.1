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

    system->updateStateEnergies();

    system->updateTransitionRates();

    system->sampleTransitionEvent();

    system->updateStateOccupation();

    system->increaseSystemTime(system->totalSumOfRates);

    if(trackingData) {
        system->eventCounts[system->lastHopIndices[0]*numOfStates + system->lastHopIndices[1]] += 1;
    }
}

void Simulator::simulateNumberOfSteps(unsigned int numOfSteps, bool trackingData) {

    for(int i = 0; i < numOfSteps; ++i) {
        singleMCStep(trackingData);
    }
}