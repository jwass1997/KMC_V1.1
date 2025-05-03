#include <algorithm>
#include <omp.h>

#include "SystemGraph.h"
#include "FEMmethods.h"

SystemGraph::SystemGraph() 
    : finiteElementSolver(nullptr)
{
    /**
     * 
     * Graph is directly ready for simulation
     * 
     */

    numOfStates = nAcceptors + nElectrodes;

    acceptorCoordinates.resize(2*nAcceptors, 0.0);
    donorCoordinates.resize(2*nDonors, 0.0);
    electrodeCoordinates.resize(2*nElectrodes, 0.0);
    distanceMatrix.resize(numOfStates*numOfStates, 0.0);
    inverseAcceptorDistances.resize(nAcceptors*nAcceptors, 0.0);
    occupationOfStates.resize(nAcceptors, 0);
    acceptorInteraction.resize(nAcceptors*nAcceptors, 0.0);
    randomEnergies.resize(nAcceptors, 0.0);
    acceptorDonorInteraction.resize(nAcceptors, 0.0);
    stateEnergies.resize(numOfStates, 0.0);    
    eventCounts.resize(numOfStates*numOfStates, 0);
    lastHopIndices.resize(2, 0);

    initializeSystemGraph();
}

SystemGraph::SystemGraph(const std::string& path) 
    : finiteElementSolver(nullptr)
{   
    initializeSystemGraph(path);
}

SystemGraph::~SystemGraph() 
{
    for(auto electrode : electrodeData) {
        delete electrode;
        electrode = nullptr;
    }
    delete finiteElementSolver;
    finiteElementSolver = nullptr;
}

const unsigned int& SystemGraph::getAcceptorNumber() const {
    return nAcceptors;
}

const unsigned int& SystemGraph::getDonorNumber() const {
    return nDonors;
}

const unsigned int& SystemGraph::getElectrodeNumber() const {
    return nElectrodes;
}

const unsigned int& SystemGraph::getNumOfStates() const {
    return numOfStates;
}

const double& SystemGraph::getDistance(unsigned int stateIndex1, unsigned int stateIndex2) const {
    return distanceMatrix[stateIndex1*numOfStates + stateIndex2];
}

const std::vector<double> SystemGraph::getAcceptorCoordinates(unsigned int acceptorIndex) const {
    return {acceptorCoordinates[acceptorIndex*2], acceptorCoordinates[acceptorIndex*2 + 1]};
}

const std::vector<double> SystemGraph::getDonorCoordinates(unsigned int donorIndex) const {
    return {donorCoordinates[donorIndex*2], donorCoordinates[donorIndex*2 + 1]};
}

const std::vector<double> SystemGraph::getElectrodeCoordinates(unsigned int electrodeIndex) const {
    return {electrodeCoordinates[electrodeIndex*2], electrodeCoordinates[electrodeIndex*2 + 1]};
}

const int& SystemGraph::getOccupationOfState(unsigned int stateIndex) const {
    return occupationOfStates[stateIndex];
}

const double& SystemGraph::getStateEnergy(unsigned int stateIndex) const {
    return stateEnergies[stateIndex];
}

const double& SystemGraph::getElectrodeVoltage(unsigned int electrodeIndex) const {
    if(electrodeIndex >= nElectrodes) {
        throw std::invalid_argument("Electrode index is out of bounds");
    }

    return electrodeData[electrodeIndex]->voltage;
}

const double& SystemGraph::getSystemTime() const {
    return systemTime;
}

const int& SystemGraph::getNumberOfEvents(int stateIndex1, int stateIndex2) const {
    return eventCounts[stateIndex1*numOfStates + stateIndex2];
}

void SystemGraph::setOccupied(unsigned int stateIndex) {
    occupationOfStates[stateIndex] = 1;
}

void SystemGraph::setUnoccupied(unsigned int stateIndex) {
    occupationOfStates[stateIndex] = 0;
}

void SystemGraph::setStateEnergy(unsigned int stateIndex, double energy) {
    stateEnergies[stateIndex] = energy;
}

void SystemGraph::setElectrodeVoltage(unsigned int electrodeIdx, double voltage) {
    if(electrodeIdx >= nElectrodes) {
        throw std::invalid_argument("Index is out of bounds");
    }

    electrodeData[electrodeIdx]->voltage = voltage;
}

void SystemGraph::initializeSystemGraph() {

    radius = radius / R;
    electrodeWidth = electrodeWidth / R;
    a = a / R;
    minHopDistance = minHopDistance / R;
    maxHopDistance = maxHopDistance / R;

    numOfStates = nAcceptors + nElectrodes;

    initializeCoordinates();
    initializeElectrodes();
    initializeContainers();
    initializeOccupiedStates();
    initializePotential();
    initializeStateEnergies();
}

std::filesystem::path SystemGraph::getConfigFilePath(const std::string& folderPath, const std::string& fileName) {
    
    return std::filesystem::path(folderPath) / fileName;
}

void SystemGraph::initializeSystemGraph(const std::string& path) {

    if(path.empty()) {
        throw std::invalid_argument("Specify a path");
    }

    initializeCoordinates(path);
    initializeElectrodes(path);

    nAcceptors = acceptorCoordinates.size() / 2;
    nDonors = donorCoordinates.size() / 2;
    nElectrodes = electrodeCoordinates.size() / 2;

    numOfStates = nAcceptors + nElectrodes;

    R = std::sqrt(M_PI*radius*radius / static_cast<double>(nAcceptors));

    radius = radius / R;
    electrodeWidth = electrodeWidth / R;
    a = a / R;
    minHopDistance = minHopDistance / R;
    maxHopDistance = maxHopDistance / R;

    for (int i = 0; i < nElectrodes; ++i) {
        electrodeCoordinates[i*2] = electrodeCoordinates[i*2] / R;
        electrodeCoordinates[i*2 + 1] = electrodeCoordinates[i*2 + 1] / R;
    }

    distanceMatrix.resize(numOfStates*numOfStates, 0.0);
    inverseAcceptorDistances.resize(nAcceptors*nAcceptors, 0.0);
    occupationOfStates.resize(nAcceptors, 0);
    randomEnergies.resize(nAcceptors, 0.0);
    acceptorDonorInteraction.resize(nAcceptors, 0.0);
    stateEnergies.resize(numOfStates, 0.0);
    acceptorInteraction.resize(nAcceptors*nAcceptors, 0.0);
    eventCounts.resize(numOfStates*numOfStates, 0);
    lastHopIndices.resize(2, 0);

    initializeContainers();
    initializeOccupiedStates();
    initializePotential();
    initializeStateEnergies();
}

void SystemGraph::initializeCoordinates() {

    for (int i = 0; i < nAcceptors; ++i) {
        double randomPhi = 2.0*M_PI*sampleFromUniformDistribution(0.0, 1.0);
        double randomR = radius*std::sqrt(sampleFromUniformDistribution(0.0, 1.0));
        acceptorCoordinates[i*2] = randomR*std::cos(randomPhi);
        acceptorCoordinates[i*2 + 1] = randomR*std::sin(randomPhi);
    }

    for (int i = 0; i < nDonors;  ++i) {
        double randomPhi = 2.0*M_PI*sampleFromUniformDistribution(0.0, 1.0);
        double randomR = radius*std::sqrt(sampleFromUniformDistribution(0.0, 1.0));
        donorCoordinates[i*2] = randomR*std::cos(randomPhi);
        donorCoordinates[i*2 + 1] = randomR*std::sin(randomPhi);
    }
}

void SystemGraph::initializeElectrodes() {

    for(int i = 0; i < nElectrodes; ++i) {
        Electrode* newElectrode = new Electrode;
        electrodeData.push_back(newElectrode);
        electrodeData[i]->angularPosition = std::get<0>(defaultCircularElectrodeConfig[i]);
        electrodeData[i]->voltage = std::get<1>(defaultCircularElectrodeConfig[i]);
    }

    for (int i = 0; i < nElectrodes; ++i) {
        double phi = (2.0*M_PI*electrodeData[i]->angularPosition) / 360.0;
        electrodeCoordinates[i*2] = radius*std::cos(phi);
        electrodeCoordinates[i*2 + 1] = radius*std::sin(phi);
    }
}

void SystemGraph::initializeCoordinates(const std::string& configuration) {

    auto acceptorConfig = getConfigFilePath(configuration, "device_circle_acceptors.txt");
    auto donorConfig = getConfigFilePath(configuration, "device_circle_donors.txt");

    std::ifstream acceptorFile(acceptorConfig);

    if(!acceptorFile.is_open()) {
        std::cerr << "No such file: " << acceptorConfig << "\n"; 
    } 
    else {   
        std::string line;

        while(getline(acceptorFile, line)) {
            if(line.empty() || line[0] == '#') {
                continue;
            }
            std::istringstream ss(line);
            double coordX, coordY;
            ss >> coordX >> coordY;
            acceptorCoordinates.push_back(coordX);
            acceptorCoordinates.push_back(coordY);
        }
        acceptorFile.close();
    }

    std::ifstream donorFile(donorConfig);

    if(!donorFile.is_open()) {
        std::cerr << "No such file: " << donorConfig << "\n"; 
    } 
    else {
        std::string line;

        while(getline(donorFile, line)) {
            if(line.empty() || line[0] == '#') {
                continue;
            }
            std::istringstream ss(line);
            double coordX, coordY;
            ss >> coordX >> coordY;
            donorCoordinates.push_back(coordX);
            donorCoordinates.push_back(coordY);
        }
        donorFile.close();
    }
}

void SystemGraph::initializeElectrodes(const std::string& configuration) {

    auto electrodeConfig = getConfigFilePath(configuration, "default_circle_electrodes.txt");
    std::ifstream electrodeFile(electrodeConfig);

    if(!electrodeFile.is_open()) {
        std::cerr << "No such file: " << configuration << "\n";
    }
    else {
        std::string line;
        while(getline(electrodeFile, line)) {
            if(line.empty() || line[0] == '#') {

            }
            else {
                std::istringstream ss(line);
                std::string angularPosition, voltage;
                ss >> angularPosition >> voltage;

                Electrode* newElectrode = new Electrode;

                newElectrode->angularPosition = std::stod(angularPosition);
                newElectrode->voltage = std::stod(voltage); 
                electrodeData.push_back(newElectrode);

                double phi = (2.0*M_PI*newElectrode->angularPosition) / 360.0;
                double x = radius*std::cos(phi);
                double y = radius*std::sin(phi);
                electrodeCoordinates.push_back(x);
                electrodeCoordinates.push_back(y);
            }
        }
        electrodeFile.close();
    }
}

void SystemGraph::initializeContainers() {

    std::vector<std::array<double,2>> allCoordinates(numOfStates);

    for (int i = 0; i < nAcceptors; ++i) {
        allCoordinates[i][0] = acceptorCoordinates[i*2];
        allCoordinates[i][1] = acceptorCoordinates[i*2 + 1];
    }

    for (int i = 0; i < nElectrodes; ++i) {
        allCoordinates[i+nAcceptors][0] = electrodeCoordinates[i*2];
        allCoordinates[i+nAcceptors][1] = electrodeCoordinates[i*2 + 1];
    }

    numOfNeighbours.resize(numOfStates);
    int totalNumOfEvents = 0;
    for (int i = 0; i < numOfStates; ++i) {
        distanceMatrix[i*numOfStates + i] = 0.0;
        for (int j = i + 1; j < numOfStates; ++j) {
            double Dx = allCoordinates[i][0] - allCoordinates[j][0];
            double Dy = allCoordinates[i][1] - allCoordinates[j][1];
            double distance = std::sqrt(Dx*Dx + Dy*Dy);
            distanceMatrix[i*numOfStates + j] = distance;
            distanceMatrix[j*numOfStates + i] = distance;
            if ((distance > minHopDistance) && (distance < maxHopDistance)) {
                totalNumOfEvents++;
                numOfNeighbours[i]+=1;
                numOfNeighbours[j]+=1;
            }
        }  
    }

    jaggedArrayLengths.resize(numOfStates+1);
    jaggedArrayLengths[0] = 0;
    for (int i = 0; i < numOfStates; ++i) {
        jaggedArrayLengths[i+1] = jaggedArrayLengths[i] + numOfNeighbours[i];
    }

    std::vector<int> writePtr(numOfStates);
    for (int i = 0; i < numOfStates; ++i) {
        writePtr[i] = jaggedArrayLengths[i];
    }

    constantTransitionRates.resize(2*totalNumOfEvents);
    dynamicalTransitionRates.resize(2*totalNumOfEvents);
    aggregatedTransitionRates.resize(2*totalNumOfEvents);
    neighbourIndices.resize(2*totalNumOfEvents);

    for (int i = 0; i < numOfStates; ++i) {
        for (int j = i+1; j < numOfStates; ++j) {
            double distance =  distanceMatrix[i*numOfStates + j];
            if (distance > minHopDistance && distance < maxHopDistance) {
                int indexIJ = writePtr[i]++;
                int indexJI = writePtr[j]++;
                neighbourIndices[indexIJ] = j;
                neighbourIndices[indexJI] = i;
                constantTransitionRates[indexIJ] = nu0*fastExp(-2.0*distance / a);
                constantTransitionRates[indexJI] = nu0*fastExp(-2.0*distance / a);
            }
        }
    }

    for (int i = 0; i < nAcceptors; ++i) {
        inverseAcceptorDistances[i*nAcceptors + i] = 0.0;
        for (int j = i+1; j < nAcceptors; ++j) {
            double inverseDistance = 1.0 / distanceMatrix[i*numOfStates + j];
            inverseAcceptorDistances[i*nAcceptors + j] = inverseDistance;
            inverseAcceptorDistances[j*nAcceptors + i] = inverseDistance;
        }
    } 
}

void SystemGraph::initializeStateEnergies() {

    std::vector<double> inverseDistances(nAcceptors, 0.0);
    // Acc-Don interaction + random energy + potential energy (for acceptors only)
    for(int i = 0; i < nAcceptors; ++i) {
        double potentialEnergy = finiteElementSolver->getPotential(
            acceptorCoordinates[i*2], 
            acceptorCoordinates[i*2 + 1]
        )*e / kbT;
        std::cout << potentialEnergy << "\n";
		double sumOfInverseDistances = 0.0;
		for(int j = 0; j <  nDonors; j++) {
			sumOfInverseDistances += 1.0 / calculateDistance(
                acceptorCoordinates[i*2], 
                donorCoordinates[j*2], 
                acceptorCoordinates[i*2 + 1], 
                donorCoordinates[j*2 + 1]
            );
		}

		inverseDistances[i] = sumOfInverseDistances;
		acceptorDonorInteraction[i] = A0*sumOfInverseDistances;

		if(addRandomEnergy) {
			std::normal_distribution<double> normalDist(0.0, energyDisorder);
			double randomEnergy = sampleFromNormalDistribution(0.0, energyDisorder);	
			randomEnergies[i] = randomEnergy;		
		}
        stateEnergies[i] = potentialEnergy + acceptorDonorInteraction[i] + randomEnergies[i];
	}
    // Potential energy (for electrodes only)
	for(int i = nAcceptors; i < nAcceptors+nElectrodes; ++i) {
		stateEnergies[i] = electrodeData[i-nAcceptors]->voltage*e / kbT;
	}
    // Acc-Acc interaction
    for (int i = 0; i < nAcceptors; ++i) {
        for (int j = 0; j < nAcceptors; ++j) {
            if (i != j) {
                acceptorInteraction[i] += (1 - occupationOfStates[j]) * inverseAcceptorDistances[i*nAcceptors + j];
            }
        }
        stateEnergies[i] += - A0*acceptorInteraction[i];
    }
}

void SystemGraph::initializeOccupiedStates() {

    if (nAcceptors <= nDonors) {
        throw std::invalid_argument("Number of acceptors can not be equal or smaller than number of donors");
    }

    std::vector<int> randomVector(nAcceptors, 0);

    for(int i = 0; i < nAcceptors; ++i) {
        randomVector[i] = i;
    }
    std::shuffle(randomVector.begin(), randomVector.end(), rng);
    for (int i = 0; i < nAcceptors - nDonors; ++i) {
        occupationOfStates[randomVector[i]] = 1;
    }
}

void SystemGraph::initializePotential() {

    finiteElementSolver = new FiniteElementeCircle(radius, 1e5);

    for (int i = 0; i < electrodeData.size(); ++i) {
        finiteElementSolver->setElectrode(
            electrodeData[i]->voltage,
            electrodeData[i]->angularPosition / 360.0 * 2.0*M_PI - 0.5*electrodeWidth,
            electrodeData[i]->angularPosition / 360.0 * 2.0*M_PI + 0.5*electrodeWidth
        );
    }

    finiteElementSolver->initRun(true);

    finiteElementSolver->run();
}

void SystemGraph::updateStateEnergies() {

    if (lastHopIndices[0] < nAcceptors && lastHopIndices[1] < nAcceptors) {
        for (int i = 0; i < nAcceptors; ++i) {
            if (i != lastHopIndices[1]) {
                acceptorInteraction[i] -= 1.0*inverseAcceptorDistances[i*nAcceptors + lastHopIndices[1]];
            }
            if (i != lastHopIndices[0]) {
                acceptorInteraction[i] += 1.0*inverseAcceptorDistances[i*nAcceptors + lastHopIndices[0]];
            }
        }
    }
    if (lastHopIndices[0] < nAcceptors && lastHopIndices[1] >= nAcceptors) {
        for (int i = 0; i < nAcceptors; ++i) {
            if (i != lastHopIndices[0]) {
                acceptorInteraction[i] += 1.0*inverseAcceptorDistances[i*nAcceptors + lastHopIndices[0]];
            }
        }
    }
    if (lastHopIndices[0] >= nAcceptors && lastHopIndices[1] < nAcceptors) {
        for (int i = 0; i < nAcceptors; ++i) {
            if (i != lastHopIndices[1]) {
                acceptorInteraction[i] -= 1.0*inverseAcceptorDistances[i*nAcceptors + lastHopIndices[1]];
            }
        }
    }

    for (int i = 0; i < nAcceptors; ++i) {
        double potentialEnergy = finiteElementSolver->getPotential(
            acceptorCoordinates[i*2],
            acceptorCoordinates[i*2 + 1]
        );
        stateEnergies[i] = potentialEnergy + acceptorDonorInteraction[i] + randomEnergies[i] - A0*acceptorInteraction[i];
    }
}

void SystemGraph::updateStateOccupation() {

    if (lastHopIndices[0] < nAcceptors && lastHopIndices[1] < nAcceptors) {
		occupationOfStates[lastHopIndices[0]] = 0;
		occupationOfStates[lastHopIndices[1]] = 1;
	}
	if (lastHopIndices[0] < nAcceptors && lastHopIndices[1] >= nAcceptors) {
		occupationOfStates[lastHopIndices[0]] = 0;
	}
	if (lastHopIndices[0] >= nAcceptors) {
		if(lastHopIndices[1] < nAcceptors) {
			occupationOfStates[lastHopIndices[1]] = 1;
		}
	}
}

void SystemGraph::updateTransitionRates() {

    double sum = 0.0;

    #pragma omp parallel for schedule(static) reduction(+:sum)

    for (int i = 0; i < numOfStates; ++i) {
		for (int k = jaggedArrayLengths[i]; k < jaggedArrayLengths[i+1]; ++k) {
            int partner = neighbourIndices[k];
			// Electrode - Acceptor
			if (i >= nAcceptors && partner < nAcceptors) {
				if(occupationOfStates[partner] == 1) {
					dynamicalTransitionRates[k] = 0.0;
				}
				else {
					double deltaE = stateEnergies[partner] - stateEnergies[i];
					if (deltaE < 0.0) {
						dynamicalTransitionRates[k] = nu0;
					} 
					else {
						dynamicalTransitionRates[k] = nu0*fastExp(-deltaE);
					} 
				}
			}
			// Acceptor - Electrode
			else if (i < nAcceptors && partner >= nAcceptors) {
				if (occupationOfStates[i] == 0) {
					dynamicalTransitionRates[k] = 0.0;
				} 
				else {
					double deltaE = stateEnergies[partner] - stateEnergies[i];
					if (deltaE < 0.0) {
						dynamicalTransitionRates[k] = nu0;
					}
					else {
						dynamicalTransitionRates[k] = nu0*fastExp(-deltaE);
					}
				}
			}
			// Acceptor - Acceptor
			else if (i < nAcceptors && partner < nAcceptors) {
				if ((occupationOfStates[i] == 1) && (occupationOfStates[partner] == 0)) {
					double deltaE = stateEnergies[partner] - stateEnergies[i] - A0 / distanceMatrix[i*numOfStates + partner];
					if (deltaE < 0.0) {
						dynamicalTransitionRates[k] = nu0;
					}
					else {
						dynamicalTransitionRates[k] = nu0*fastExp(-deltaE);
					} 
				}
				else {
					dynamicalTransitionRates[k] = 0.0;
				}					
			}

            aggregatedTransitionRates[k] = dynamicalTransitionRates[k]*constantTransitionRates[k];
            sum += aggregatedTransitionRates[k];
		}
	}

    totalSumOfRates = sum;
}

void SystemGraph::sampleTransitionEvent() {

    double r = sampleFromUniformDistribution(0.0, totalSumOfRates);
    cumulativeSumOfRates = 0.0;
    for (int i = 0; i < numOfStates; ++i) {
        int L = jaggedArrayLengths[i];
        int R = jaggedArrayLengths[i+1];
        for (int k = L; k < R; ++k) {
            double rate = aggregatedTransitionRates[k];
            cumulativeSumOfRates+=rate;
            if (cumulativeSumOfRates >= r) {
                int j = neighbourIndices[k];
                lastHopIndices[0] = i;
                lastHopIndices[1] = j;
                return;
            }
        }
    }
}

void SystemGraph::increaseEventCount() {
    eventCounts[lastHopIndices[0]*numOfStates + lastHopIndices[1]] += 1;
}

void SystemGraph::increaseSystemTime(double totalSumOfTransitionRates) {

    double uniformSample = sampleFromUniformDistribution(0.0, 1.0) + 1e-12;
	double timeStep = -(std::log(uniformSample) / totalSumOfTransitionRates);

    systemTime += timeStep;
}

void SystemGraph::resetEventCounts() {
    std::fill(eventCounts.begin(), eventCounts.end(), 0);
}

void SystemGraph::resetPotential() {
    // Subtracts current potential energy
    for (int i = 0; i < nAcceptors; ++i) {
        stateEnergies[i] -= finiteElementSolver->getPotential(acceptorCoordinates[i*2], acceptorCoordinates[i*2 + 1]);
    }
}

void SystemGraph::updatePotential() {
    // Adds current potential energy
    for (int i = 0; i < nAcceptors; ++i) {
        stateEnergies[i] += finiteElementSolver->getPotential(acceptorCoordinates[i*2], acceptorCoordinates[i*2 + 1]);
    }
}

void SystemGraph::updateVoltages(std::vector<double>& voltages) {
    /**
     * 
     * Whole process of updating an electrode and updating constant state energies.
     * After calling this function the device is ready to be simulated again.
     * 
     */

    if(finiteElementSolver == nullptr) {
        std::cerr << "Solver not initialized, cant update potential" << "\n";
    }
    else if(voltages.size() > nElectrodes) {
        throw std::invalid_argument("Too many voltages");
    }

    resetPotential();

    for (int i = 0; i < voltages.size(); ++i) {
        finiteElementSolver->updateElectrodeVoltage(i, voltages[i]);
    }

    finiteElementSolver->run();
        
    updatePotential();

    for (int i = 0; i < nElectrodes; ++i) {
        stateEnergies[i + nAcceptors] = voltages[i]*e / kbT;
    }
}