#pragma once

#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <tuple>
#include <random>
#include <string>
#include <filesystem>
#include <sstream>
#include <fstream>

#include "utils.h"

class CircularFEMSolver;

class SystemGraph {
    public:
        SystemGraph();

        SystemGraph(const std::string& path);

        ~SystemGraph();

        std::filesystem::path getConfigFilePath(const std::string& folderPath, const std::string& fileName);

        void initializeSystemGraph();

        void initializeSystemGraph(const std::string& path);
   
        void initializeCoordinates();

        void initializeElectrodes();

        void initializeCoordinates(const std::string& configuration);

        void initializeElectrodes(const std::string& configuration);

        void initializeContainers();

        void initializeStateEnergies();

        void initializeOccupiedStates();

        void initializePotential();
        
        const unsigned int& getAcceptorNumber() const;

        const unsigned int& getDonorNumber() const;

        const unsigned int& getElectrodeNumber() const;

        const unsigned int& getNumOfStates() const;

        const double& getDistance(unsigned int stateIndex1, unsigned int stateIndex2) const;

        const std::vector<double> getAcceptorCoordinates(unsigned int acceptorIndex) const;

        const std::vector<double> getDonorCoordinates(unsigned int donorIndex) const;

        const std::vector<double> getElectrodeCoordinates(unsigned int electrodeIndex) const;

        const int& getOccupationOfState(unsigned int stateIndex) const;

        const double& getStateEnergy(unsigned int stateIndex) const;

        const double& getElectrodeVoltage(unsigned int electrodeIndex) const;

        const double& getSystemTime() const;

        const int& getNumberOfEvents(int stateIndex1, int stateIndex2) const;

        void setOccupied(unsigned int stateIndex);

        void setUnoccupied(unsigned int stateIndex);

        void setStateEnergy(unsigned int stateIndex, double energy);

        void setElectrodeVoltage(unsigned int electrodeIdx, double voltage);

        void updateStateEnergies();

        void updateStateOccupation();

        void updateTransitionRates();

        int selectEvent(std::vector<double> cumulativeRateCatalog);

        void sampleTransitionEvent();

        void updateElectrodeVoltage(int electrodeIndex, double voltage);

        void multiElectrodeUpdate(std::vector<int> electrodeIndices, std::vector<double> voltages);

        void increaseEventCount();

        void increaseSystemTime(double totatlSumOfTransitionRates);

        void resetEventCounts();

        void resetPotential();

        void updatePotential();

        friend class Simulator;

    private:
        bool dimensionlessDistances;

        unsigned int nAcceptors = 200;
        unsigned int nDonors = 3;
        unsigned int nElectrodes = 8;
    
        unsigned int numOfStates;

        double radius = 100.0;

        double nu0 = 1.0;
        double a = 20.0;
        double T = 77.0;

        double energyDisorder = 0.05*e / kbT;
        double R = std::sqrt(M_PI*radius*radius / static_cast<double>(nAcceptors));
        double A0 = (e*e) / (4.0*kbT*PI*eps0*epsr*1e-9*R);

        bool addRandomEnergy = true;

        double minHopDistance = 3.0;

        double maxHopDistance = 60.0;

        double systemTime = 0.0;

        double totalSumOfRates;

        double cumulativeSumOfRates;

        std::vector<double> acceptorCoordinates;

        std::vector<double> donorCoordinates;

        std::vector<double> electrodeCoordinates;

        std::vector<double> distanceMatrix;

        std::vector<double> inverseAcceptorDistances;

        std::vector<Electrode* > electrodeData;

        std::vector<int> eventCounts;

        std::vector<int> occupationOfStates;

        std::vector<double> constantStateEnergies;

        std::vector<double> stateEnergies;

        std::vector<double> acceptorInteraction;

        std::vector<double> constantTransitionRates;

        std::vector<double> dynamicalTransitionRates;

        std::vector<double> aggregatedTransitionRates;

        std::vector<int> numOfNeighbours;
        
        std::vector<int> jaggedArrayLengths;

        std::vector<int> neighbourIndices;

        std::vector<int> lastHopIndices;

        std::vector<std::tuple<double, double>> defaultCircularElectrodeConfig = {
            {0.0, -1.0},
            {45.0, 1.0},
            {90.0, -1.0},
            {135.0, 1.0},
            {180.0, -1.0},
            {225.0, 1.0},
            {270.0, -1.0},
            {315.0, 1.0},
        };

        CircularFEMSolver* finiteElementSolver;

        friend class Simulator;
};