#pragma once

#include <filesystem>

#include "utils.h"

class SystemGraph;

class Simulator {
    public:
        Simulator();

        Simulator(const std::string& path);

        Simulator(const std::string& path, int burnInHops, int numOfHops);

        ~Simulator();

        void singleMCStep(bool trackingData);

        void simulateNumberOfSteps(unsigned int numOfSteps, bool trackingData);

        SystemGraph* system;

    private:
        unsigned int burnInHops = 1e4;

        unsigned int numOfHops = 1e5;

        double simulationTime = 0.0;
};