#include "SystemGraph.h"
#include "CircularFEMSolver.h"

SystemGraph::SystemGraph() 
    : acceptorCoordinates(nAcceptors, std::vector<double>(2, 0.0))
    , donorCoordinates(nDonors, std::vector<double>(2, 0.0))
    , electrodeCoordinates(nElectrodes, std::vector<double>(2, 0.0))
    , distanceMatrix(nAcceptors+nElectrodes, std::vector<double>(nAcceptors+nElectrodes, 0.0))
    , occupationOfStates(nAcceptors, 0)
    , constantStateEnergies(nAcceptors+nElectrodes, 0.0)
    , stateEnergies(nAcceptors+nElectrodes, 0.0)
    , constantTransitionRates(nAcceptors+nElectrodes, std::vector<double>())
    , dynamicalTransitionRates(nAcceptors+nElectrodes, std::vector<double>())
    , listOfHoppingPartners(nAcceptors+nElectrodes, std::vector<int>())
    , eventCounts(nAcceptors+nElectrodes, std::vector<int>(nAcceptors+nElectrodes, 0))
    , finiteElementSolver(nullptr)
{
    /**
     * 
     * Graph is directly ready for simulation
     * 
     */

    initializeSystemGraph();
}

SystemGraph::SystemGraph(const std::string& path) 
    : finiteElementSolver(nullptr)
{   
    /**
     * 
     * Initialization of the system from a configuration 
     * 
     */
    
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
    return distanceMatrix[stateIndex1][stateIndex2];
}

const std::vector<double>& SystemGraph::getAcceptorCoordinates(unsigned int acceptorIndex) const {
    return acceptorCoordinates[acceptorIndex];
}

const std::vector<double>& SystemGraph::getDonorCoordinates(unsigned int donorIndex) const {
    return donorCoordinates[donorIndex];
}

const std::vector<double>& SystemGraph::getElectrodeCoordinates(unsigned int electrodeIndex) const {
    return electrodeCoordinates[electrodeIndex];
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
    return eventCounts[stateIndex1][stateIndex2];
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
    a = a / R;
    minHopDistance = minHopDistance / R;
    maxHopDistance = maxHopDistance / R;

    numOfStates = nAcceptors + nElectrodes;

    initializeCoordinates();
    initializeElectrodes();
    initializeContainers();
    initializeOccupiedStates();
    initializePotential();
    initializeConstantStateEnergies();
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

    nAcceptors = acceptorCoordinates.size();
    nDonors = donorCoordinates.size();
    nElectrodes = electrodeCoordinates.size();

    numOfStates = nAcceptors + nElectrodes;

    R = std::sqrt(M_PI*radius*radius / static_cast<double>(nAcceptors));

    radius = radius / R;
    a = a / R;
    minHopDistance = minHopDistance / R;
    maxHopDistance = maxHopDistance / R;

    for (auto& coordinates : electrodeCoordinates) {
        coordinates[0] = coordinates[0] / R;
        coordinates[1] = coordinates[1] / R;
    }

    distanceMatrix.resize(nAcceptors+nElectrodes, std::vector<double>(nAcceptors+nElectrodes, 0.0));
    occupationOfStates.resize(nAcceptors, 0);
    constantStateEnergies.resize(nAcceptors+nElectrodes, 0.0);
    stateEnergies.resize(nAcceptors+nElectrodes, 0.0);
    constantTransitionRates.resize(nAcceptors+nElectrodes, std::vector<double>());
    dynamicalTransitionRates.resize(nAcceptors+nElectrodes, std::vector<double>());
    listOfHoppingPartners.resize(nAcceptors+nElectrodes, std::vector<int>());
    eventCounts.resize(nAcceptors+nElectrodes, std::vector<int>(nAcceptors+nElectrodes, 0));
    initializeContainers();
    initializeOccupiedStates();
    initializePotential();
    initializeConstantStateEnergies();
}

void SystemGraph::initializeCoordinates() {

    /**
     * 
     * Default initialization
     * 
     */

    for (auto& coordinates : acceptorCoordinates) {
        double randomPhi = 2.0*M_PI*sampleFromUniformDistribution(0.0, 1.0);
        double randomR = radius*std::sqrt(sampleFromUniformDistribution(0.0, 1.0));

        coordinates[0] = randomR*std::cos(randomPhi);
        coordinates[1] = randomR*std::sin(randomPhi);
    }

    for (auto& coordinates : donorCoordinates) {
        double randomPhi = 2.0*M_PI*sampleFromUniformDistribution(0.0, 1.0);
        double randomR = radius*std::sqrt(sampleFromUniformDistribution(0.0, 1.0));

        coordinates[0] = randomR*std::cos(randomPhi);
        coordinates[1] = randomR*std::sin(randomPhi);
    }
}

void SystemGraph::initializeElectrodes() {

    /**
     * 
     * Default initialization
     * 
     */

    for(int i = 0; i < nElectrodes; ++i) {
        Electrode* newElectrode = new Electrode;
        electrodeData.push_back(newElectrode);
        electrodeData[i]->angularPosition = std::get<0>(defaultCircularElectrodeConfig[i]);
        electrodeData[i]->voltage = std::get<1>(defaultCircularElectrodeConfig[i]);
    }

    for (int i = 0; i < nElectrodes; ++i) {
        double phi = (2.0*M_PI*electrodeData[i]->angularPosition) / 360.0;
        electrodeCoordinates[i][0] = radius*std::cos(phi);
        electrodeCoordinates[i][1] = radius*std::sin(phi);
    }
}

void SystemGraph::initializeCoordinates(const std::string& configuration) {

    /**
     * 
     * Initialize from config
     * 
     */

    auto acceptorConfig = getConfigFilePath(configuration, "default_circle_acceptors.txt");
    auto donorConfig = getConfigFilePath(configuration, "default_circle_donors.txt");

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
            acceptorCoordinates.push_back({coordX, coordY});
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
            donorCoordinates.push_back({coordX, coordY});
        }
        donorFile.close();
    }
}

void SystemGraph::initializeElectrodes(const std::string& configuration) {

    /**
     * 
     * Initialize from config
     * 
     */

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
                electrodeCoordinates.push_back({x, y});
            }
        }
        electrodeFile.close();
    }
}

void SystemGraph::initializeContainers() {

    std::vector<std::vector<double>> allCoordinates;

    for (const auto& nodeCoordinates : acceptorCoordinates) {
        allCoordinates.push_back(nodeCoordinates);
    }

    for (const auto& nodeCoordinates : electrodeCoordinates) {
        allCoordinates.push_back(nodeCoordinates);
    }

    int totalNumOfEvents = 0;
    for (int i = 0; i < nAcceptors + nElectrodes; ++i) {
        distanceMatrix[i][i] = 0.0;
        for (int j = i + 1; j < nAcceptors + nElectrodes; ++j) {
            double Dx = allCoordinates[i][0] - allCoordinates[j][0];
            double Dy = allCoordinates[i][1] - allCoordinates[j][1];
            double distance = std::sqrt(Dx*Dx + Dy*Dy);
            distanceMatrix[i][j] = distance;
            distanceMatrix[j][i] = distance;
            if ((distance > minHopDistance) && (distance < maxHopDistance)) {
                listOfHoppingPartners[i].push_back(j);
                listOfHoppingPartners[j].push_back(i);
                double constRate = nu0*std::exp(-2.0*distance / a);
                constantTransitionRates[i].push_back(constRate);
                constantTransitionRates[j].push_back(constRate);
                totalNumOfEvents++;
            }
        }
        dynamicalTransitionRates[i].resize(constantTransitionRates[i].size());      
    }
    //std::cout<<totalNumOfEvents<<"\n";
}

void SystemGraph::initializeConstantStateEnergies() {

    /**
     * 
     * Repulsion and potential energy
     * 
     */

    std::vector<double> inverseDistances(nAcceptors, 0.0);

    for(int i = 0; i < nAcceptors; ++i) {
        double energy = finiteElementSolver->calculatePotential(
            acceptorCoordinates[i][0], 
            acceptorCoordinates[i][1]
        )*e / kbT;
		double sumOfInverseDistances = 0.0;
		for(int j = 0; j <  nDonors; j++) {
			sumOfInverseDistances += 1.0 / calculateDistance(
                acceptorCoordinates[i][0], 
                donorCoordinates[j][0], 
                acceptorCoordinates[i][1], 
                donorCoordinates[j][1]
            );
		}

		inverseDistances[i] = sumOfInverseDistances;
		constantStateEnergies[i] = A0*sumOfInverseDistances;

		constantStateEnergies[i] += energy;
		if(addRandomEnergy) {
			std::normal_distribution<double> normalDist(0.0, energyDisorder);
			double randomEnergy = sampleFromNormalDistribution(0.0, energyDisorder);	
			constantStateEnergies[i] += randomEnergy;		
		}
	}

	for(int i = nAcceptors; i < nAcceptors+nElectrodes; ++i) {
		constantStateEnergies[i] = electrodeData[i-nAcceptors]->voltage*e / kbT;
	}
}

void SystemGraph::initializeOccupiedStates() {

    std::vector<int> randomVector(nAcceptors, 0);

    for(int i = 0; i < nAcceptors; ++i) {
        randomVector[i] = i;
    }
    std::shuffle(randomVector.begin(), randomVector.end(), rng);
    for(int i = 0; i < nAcceptors - nDonors; ++i) {
        occupationOfStates[randomVector[i]] = 1;
    }
}

void SystemGraph::initializePotential() {

    finiteElementSolver = new CircularFEMSolver(radius, 200, 200, electrodeData);
    finiteElementSolver->solvePoissonEquation();
}

void SystemGraph::updateStateEnergy(unsigned int indexOfStateToUpdate) {

    if(indexOfStateToUpdate >= nAcceptors+nElectrodes) {
        throw std::invalid_argument("State index out of bounds !");
    }

    if(indexOfStateToUpdate < nAcceptors) {
        double acceptorInteraction = 0.0;
        for(int j = 0; j < nAcceptors; ++j) {
            if(j != indexOfStateToUpdate) {
                acceptorInteraction += (1 - occupationOfStates[j]) / getDistance(indexOfStateToUpdate, j);
            }
        }
        stateEnergies[indexOfStateToUpdate] = constantStateEnergies[indexOfStateToUpdate];
        stateEnergies[indexOfStateToUpdate] -= A0*acceptorInteraction;
    }

    if(indexOfStateToUpdate >= nAcceptors) {
        stateEnergies[indexOfStateToUpdate] = electrodeData[indexOfStateToUpdate-nAcceptors]->voltage*e / kbT;
    }
}

void SystemGraph::updateStateOccupation(int fromStateIndex, int toStateIndex) {

    if(fromStateIndex < nAcceptors && toStateIndex < nAcceptors) {
		occupationOfStates[fromStateIndex] = 0;
		occupationOfStates[toStateIndex] = 1;
	}
	if(fromStateIndex < nAcceptors && toStateIndex >= nAcceptors) {
		occupationOfStates[fromStateIndex] = 0;
	}
	if(fromStateIndex >= nAcceptors) {
		if(toStateIndex < nAcceptors) {
			occupationOfStates[toStateIndex] = 1;
		}
	}
}

void SystemGraph::updateTransitionRates() {

    for (int i = 0; i < listOfHoppingPartners.size(); ++i) {
		for (int j = 0; j < listOfHoppingPartners[i].size(); ++j) {
            int partner = listOfHoppingPartners[i][j];
            // Electrode - Electrode
            if (i >= nAcceptors && partner >= nAcceptors) {
                dynamicalTransitionRates[i][j] = 0.0;
            }
			// Electrode - Acceptor
			else if (i >= nAcceptors && partner < nAcceptors) {
				if(getOccupationOfState(partner) == 1) {
					dynamicalTransitionRates[i][j] = 0.0;
				}
				else {
					double deltaE = getStateEnergy(partner) - getStateEnergy(i);
					if (deltaE < 0.0) {
						dynamicalTransitionRates[i][j] = nu0;
					} 
					else {
						dynamicalTransitionRates[i][j] = nu0*std::exp(-deltaE);
					} 
				}
			}
			// Acceptor - Electrode
			else if (i < nAcceptors && partner >= nAcceptors) {
				if (getOccupationOfState(i) == 0) {
					dynamicalTransitionRates[i][j] = 0.0;
				} 
				else {
					double deltaE = getStateEnergy(partner) - getStateEnergy(i);
					if (deltaE < 0.0) {
						dynamicalTransitionRates[i][j] = nu0;
					}
					else {
						dynamicalTransitionRates[i][j] = nu0*std::exp(-deltaE);
					}
				}
			}
			// Acceptor - Acceptor
			else if (i < nAcceptors && partner < nAcceptors) {
				if ((getOccupationOfState(i) == 1) && (getOccupationOfState(partner) == 0)) {
					double deltaE = getStateEnergy(partner) - getStateEnergy(i) - A0 / getDistance(i, partner);
					if (deltaE < 0.0) {
						dynamicalTransitionRates[i][j] = nu0;
					}
					else {
						dynamicalTransitionRates[i][j] = nu0*std::exp(-deltaE);
					} 
				}
				else {
					dynamicalTransitionRates[i][j] = 0.0;
				}					
			}
		}
	}
}

int SystemGraph::selectEvent(std::vector<double> cumulativeRateCatalog) {

    double _event = sampleFromUniformDistribution(0.0, 1.0);
    int L = 0;
    int R = cumulativeRateCatalog.size() - 1;
    while(L <= R) {
        int M = static_cast<int>(L + std::floor(double(R-L) / 2.0));
        if(cumulativeRateCatalog[M] < _event) {
            L = M + 1;
        }
        else if(cumulativeRateCatalog[M] > _event) {
            R = M - 1;
        }
    }
 
    return L;
}

std::tuple<int, int, double> SystemGraph::sampleTransitionEvent() {

    std::vector<std::tuple<int, int, double>> validEvents;
    double totalRate = 0.0;

    for (int i = 0; i < numOfStates; ++i) {
        for (int j = 0; j < listOfHoppingPartners[i].size(); ++j) {
            int hoppingPartner = listOfHoppingPartners[i][j];
            double rate = dynamicalTransitionRates[i][j] * constantTransitionRates[i][j];
            if (rate > 0.0) {
                validEvents.push_back(std::make_tuple(i, hoppingPartner, rate));
                totalRate += rate;
            }
        }
    }

    double r = sampleFromUniformDistribution(0.0, totalRate);
    double cumulativeTransitionRates = 0.0;

    int selectedI = -1;
    int selectedJ = -1;

    for (const auto& event : validEvents) {
        cumulativeTransitionRates += std::get<2>(event);
        if (cumulativeTransitionRates >= r) {
            selectedI = std::get<0>(event);
            selectedJ = std::get<1>(event);    
            break;          
        }
    }

    if (selectedI == -1 || selectedJ == -1) {
        std::cerr << "No valid event indices found" << "\n";
        exit(1);
    }

    auto outputTuple = std::make_tuple(selectedI, selectedJ, totalRate);

    return outputTuple;
}

void SystemGraph::increaseEventCount(unsigned int stateIndex1, unsigned int stateIndex2) {
    eventCounts[stateIndex1][stateIndex2] += 1;
}

void SystemGraph::increaseSystemTime(double totalSumOfTransitionRates) {

    double uniformSample = sampleFromUniformDistribution(0.0, 1.0) + 1e-12;
	double timeStep = -(std::log(uniformSample) / totalSumOfTransitionRates);

    systemTime += timeStep;
}

void SystemGraph::resetEventCounts() {
    for(auto& row : eventCounts) {
        std::fill(row.begin(), row.end(), 0);
    }
}

void SystemGraph::resetPotential() {

    for (int i = 0; i < nAcceptors; ++i) {
        constantStateEnergies[i] -= finiteElementSolver->calculatePotential(acceptorCoordinates[i][0], acceptorCoordinates[i][1]);
    }
}

void SystemGraph::updatePotential() {

    for (int i = 0; i < nAcceptors; ++i) {
        constantStateEnergies[i] += finiteElementSolver->calculatePotential(acceptorCoordinates[i][0], acceptorCoordinates[i][0]);
    }
}

void SystemGraph::updateElectrodeVoltage(int electrodeIndex, double voltage) {
    /**
     * 
     * Whole process of updating an electrode and updating constant state energies.
     * After calling this function the device is ready to be simulated again.
     * 
     */

    if(finiteElementSolver == nullptr) {
        std::cerr << "Solver not initialized, cant update potential" << "\n";
    }
    else if(electrodeIndex >= nElectrodes) {
        throw std::invalid_argument("Electrode index out of bounds");
    }

    setElectrodeVoltage(electrodeIndex, voltage);
    for(int i = nAcceptors; i < nAcceptors+nElectrodes; ++i) {
        constantStateEnergies[i] = electrodeData[i-nAcceptors]->voltage*e / kbT;
    }

    double angularPosition = electrodeData[electrodeIndex]->angularPosition;

    resetPotential();

    finiteElementSolver->updateElectrode(angularPosition, electrodeIndex, electrodeData[electrodeIndex]->voltage);
    finiteElementSolver->solvePoissonEquation();
    
    updatePotential();
}

void SystemGraph::multiElectrodeUpdate(std::vector<int> electrodeIndices, std::vector<double> newVoltages) {

    if (finiteElementSolver == nullptr) {
        std::cerr << "Solver not initialized, can not update potential" << "\n";
    }

    for (const int& electrodeIndex : electrodeIndices) {
        if (electrodeIndex >= nElectrodes || electrodeIndex < 0) {
            throw std::invalid_argument("Electrode index is out of bounds");
        }
    }

    if (electrodeIndices.size() != newVoltages.size()) {
        throw std::invalid_argument("Number of electrode indices do not match voltage indices");
    }
    
    for (int i = 0; i < electrodeIndices.size(); ++i) {
        setElectrodeVoltage(electrodeIndices[i], newVoltages[i]);
        constantStateEnergies[nAcceptors + electrodeIndices[i]] = newVoltages[i]*e / kbT;
    }

    resetPotential();

    finiteElementSolver->updateMultipleElectrodes(electrodeIndices, newVoltages);
    finiteElementSolver->solvePoissonEquation();

    updatePotential();
}