#include <random>
#include <omp.h>

#include "utils.h"
#include "SystemGraph.h"
#include "Simulator.h"

std::random_device randomDevice;
std::mt19937 rng((randomDevice()));

double sampleFromUniformDistribution(double min, double max) {

    /*
    
        Draw uniformly between min-max
    
    */

    std::uniform_real_distribution<double> uniformDist(min, max);
    double sample = uniformDist(rng);

    return sample;
}

double sampleFromNormalDistribution(double mean, double standardDeviation) {

    /*

        Returns sample from a normal distribution with mean and standardDeviation

    */

    std::normal_distribution<double> normalDistribution(mean, standardDeviation);

    double sample = normalDistribution(rng);

    return sample;
}

double calculateDistance(
    double coordinateX1, 
    double coordinateX2, 
    double coordinateY1, 
    double coordinateY2) {

    double Dx = coordinateX2 - coordinateX1;
    double Dy = coordinateY2 - coordinateY1;
    double distance = std::sqrt(Dx*Dx + Dy*Dy);
    
    return distance;
}

void createDirectoryFromStringPath(const std::string& path, const std::string& folderName) {

    /**
     * 
     * Checks if directory already exists and creates it if not. The current path is used per default.
     * 
     */

    if(folderName.empty()) {
        std::cerr << "Must specify a folder name " << "\n";
        return;
    }

    std::filesystem::path directoryPath = path.empty() ? std::filesystem::current_path() : std::filesystem::path(path);
    std::filesystem::path newFolder = std::filesystem::path(directoryPath)/folderName;

    try {
        if(std::filesystem::create_directory(newFolder)) {
            std::cout << "Folder has been created: " << newFolder << "\n";
        }
        else {
            std::cout << "Folder already exists: " << newFolder << "\n";
        }
    }
    catch(const std::filesystem::filesystem_error& e) {
        std::cerr << "Error creating new folder: " << e.what() << "\n";
    }
}

void recordDevice(
    const std::string& ID, 
    int equilibriumSteps, 
    int numOfSteps, 
    const std::string& deviceConfigs, 
    const std::string& saveFolderPath) {

    /**
     * 
     * Folder to save the batch needs to be specified.
     * 
     */

    if(saveFolderPath.empty()) {
        throw std::invalid_argument("No save folder specified !");
    }

    Simulator simulator(deviceConfigs);

    int nAcceptors = simulator.system->getAcceptorNumber();
    int nElectrodes = simulator.system->getElectrodeNumber();
    int nDonors = simulator.system->getDonorNumber();

    std::vector<double> flattenedAcceptorCoordinates;
    std::vector<double> flattenedDonorCoordinates;
    std::vector<double> flattenedElectrodeCoordinates;    
    std::vector<int> flattenedEventCounts;

    std::vector<size_t> shapeFlattenedAcceptorCoordinates = {static_cast<size_t>(nAcceptors), 2};
    std::vector<size_t> shapeFlattenedDonorCoordinates = {static_cast<size_t>(nDonors), 2};
    std::vector<size_t> shapeFlattenedElectrodeCoordinates = {static_cast<size_t>(nElectrodes), 2};
    std::vector<size_t> shapeFlattenedEventCounts = {static_cast<size_t>(nAcceptors+nElectrodes), static_cast<size_t>(nAcceptors+nElectrodes)};

    std::string deviceName = saveFolderPath + "/device" + ID + ".npz";
    cnpy::npz_save(deviceName, "ID", &ID, {1}, "w"); 

    simulator.simulateNumberOfSteps(equilibriumSteps, false);
    simulator.simulateNumberOfSteps(numOfSteps, true);

    for(int i = 0; i < nAcceptors; ++i) {
        std::vector<double> coordinates = simulator.system->getAcceptorCoordinates(i);
        flattenedAcceptorCoordinates.push_back(coordinates[0]);
        flattenedAcceptorCoordinates.push_back(coordinates[1]);
    }
    for(int i = 0; i < nDonors; ++i) {
        std::vector<double> coordinates = simulator.system->getDonorCoordinates(i);
        flattenedDonorCoordinates.push_back(coordinates[0]);
        flattenedDonorCoordinates.push_back(coordinates[1]);
    }
    for(int i = 0; i < nElectrodes; ++i) {
        std::vector<double> coordinates = simulator.system->getElectrodeCoordinates(i);
        flattenedElectrodeCoordinates.push_back(coordinates[0]);
        flattenedElectrodeCoordinates.push_back(coordinates[1]);
    }

    for(int j = 0; j < nAcceptors+nElectrodes; ++j) {
        for(int i = 0; i <nAcceptors+nElectrodes; ++i) {
            flattenedEventCounts.push_back(simulator.system->getNumberOfEvents(j, i));
        }
    }

    double total_time = simulator.system->getSystemTime();

    cnpy::npz_save(deviceName, "acceptor_coordinates", flattenedAcceptorCoordinates.data(), shapeFlattenedAcceptorCoordinates, "a");
    cnpy::npz_save(deviceName, "donor_coordinates", flattenedDonorCoordinates.data(), shapeFlattenedDonorCoordinates, "a");
    cnpy::npz_save(deviceName, "electrode_coordinates", flattenedElectrodeCoordinates.data(), shapeFlattenedElectrodeCoordinates, "a");
    cnpy::npz_save(deviceName, "event_counts", flattenedEventCounts.data(), shapeFlattenedEventCounts, "a");
    cnpy::npz_save(deviceName, "device_time", &total_time, {1}, "a");
}

void IVCurve(
    double maxVoltage,
    double minVoltage,
    int numOfSimulations,
    int numOfVoltageSettings,
    int scanElectrodeIndex,
    int equilibriumSteps,
    int simulationSteps,
    int numOfIntervals,
    int ID,
    const std::string& defaultConfig,
    const std::string& saveFolderPath) {
        
        if (saveFolderPath.empty()) {
            throw std::invalid_argument("No save folder specified");
        }

        std::string dataFileName = saveFolderPath + "/IVcharacteristic" + std::to_string(ID) + ".npz";
        cnpy::npz_save(dataFileName, "ID", &ID, {1}, "w");

        std::vector<double> inputs(numOfVoltageSettings*8, 0.0);
        std::vector<double> outputs(numOfVoltageSettings, 0.0);
        
        #pragma omp parallel
        {      
            std::mt19937 threadRng(std::random_device{}() + omp_get_thread_num());
            std::uniform_real_distribution<double> uni(minVoltage, maxVoltage);

            #pragma omp for
            for (int setting = 0; setting < numOfVoltageSettings; ++setting) {
                std::vector<double> voltageSetting(8, 0.0);
                for (int i = 0; i < 8; ++i) {
                    voltageSetting[i] = uni(threadRng);
                }
                double averageOutputCurrent = 0.0;
                for (int simCount = 0; simCount < numOfSimulations; ++simCount) {
                    double current = currentFromVoltageCombination(
                        voltageSetting,
                        scanElectrodeIndex,
                        equilibriumSteps,
                        simulationSteps,
                        numOfIntervals,
                        defaultConfig
                    );

                    averageOutputCurrent += current / static_cast<double>(numOfSimulations);
                }
                outputs[setting] = averageOutputCurrent;
                for (int i = 0; i < 8; ++i) {
                    inputs[setting*8 + i] = voltageSetting[i];
                }
            }
            
        }
        std::vector<size_t> inputShape = {static_cast<size_t>(numOfVoltageSettings), 8};
        std::vector<size_t> outputShape = {static_cast<size_t>(numOfVoltageSettings)};
        cnpy::npz_save(dataFileName, "inputs", inputs.data(), inputShape, "a");
        cnpy::npz_save(dataFileName, "outputs", outputs.data(), outputShape, "a");
}

double currentFromVoltageCombination(
    std::vector<double> voltageSetting,
    int scanElectrodeIndex,
    int equilibriumSteps,
    int simulationSteps,
    int numOfIntervals,
    const std::string& defaultConfig) {
        
        Simulator simulator(defaultConfig);
        int nAcceptors = simulator.system->getAcceptorNumber();
        int numOfStates = simulator.system->getNumOfStates();
        std::vector<int> electrodeIndices = {0, 1, 2, 3, 4, 5, 6, 7};
        simulator.system->multiElectrodeUpdate(electrodeIndices, voltageSetting); 
        
        simulator.simulateNumberOfSteps(equilibriumSteps, false);

        double averageCurrent = 0.0;
        int intervalSteps = simulationSteps / numOfIntervals;
        int intervalCounter = 0;
        while(intervalCounter < numOfIntervals) {
            double startClock = simulator.system->getSystemTime();
            simulator.simulateNumberOfSteps(intervalSteps, true);
            double endClock = simulator.system->getSystemTime();

            double elapsedTime = endClock - startClock;
            int inCounts = 0;
            int outCounts = 0;
            for (int i = 0; i < numOfStates; ++i) {
                outCounts += simulator.system->getNumberOfEvents(nAcceptors+scanElectrodeIndex, i);
                inCounts += simulator.system->getNumberOfEvents(i, nAcceptors+scanElectrodeIndex); 
            }
            averageCurrent += static_cast<double>(inCounts-outCounts) / (elapsedTime*static_cast<double>(numOfIntervals));
            
            simulator.system->resetEventCounts();
            intervalCounter++;
        }

        return averageCurrent;
}

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
    const std::string& saveFolderPath, 
    int batchID) {

    if(saveFolderPath.empty()) {
        throw std::invalid_argument("No save folder specified !");
    }

    std::vector<double> inputs(batchSize*inputElectrodes.size(), 0.0);
    std::vector<double> outputs(batchSize*outputElectrodes.size(), 0.0);   
    std::vector<size_t> shapeInputs = {static_cast<size_t>(batchSize), inputElectrodes.size()}; 
    std::vector<size_t> shapeOutputs = {static_cast<size_t>(batchSize), outputElectrodes.size()}; 

    std::vector<int> systemElectrodes;
    systemElectrodes.insert(systemElectrodes.end(), inputElectrodes.begin(), inputElectrodes.end());
    systemElectrodes.insert(systemElectrodes.end(), outputElectrodes.begin(), outputElectrodes.end());

    int numOfInputElectrodes = inputElectrodes.size();
    int numOfOutputElectrodes = outputElectrodes.size();

    std::vector<double> voltages(numOfInputElectrodes+numOfOutputElectrodes);

    #pragma omp parallel 
    {   
        std::mt19937 rng(std::random_device{}() + omp_get_thread_num());
        std::uniform_real_distribution<double> uni(minVoltage, maxVoltage);
        #pragma omp for schedule(dynamic)
        for (int batch = 0; batch < batchSize; ++batch) {
            Simulator simulator(defaultConfigs);
            int nAcceptors = simulator.system->getAcceptorNumber();
            int numOfStates = simulator.system->getNumOfStates();

        
            for (int i = 0; i < numOfInputElectrodes; ++i) {
                voltages[i] = uni(rng);
                inputs[batch*numOfInputElectrodes + i] = voltages[i];
            }
            for (int i = numOfInputElectrodes; i < numOfInputElectrodes+numOfOutputElectrodes; ++i) {
                voltages[i] = -1.0;
            }
            simulator.system->multiElectrodeUpdate(systemElectrodes, voltages);
            simulator.simulateNumberOfSteps(equilibriumSteps, false);

            int intervalSteps = simulationSteps / numOfIntervals;
            int intervalCounter = 0;
            while (intervalCounter < numOfIntervals) {
                double startClock = simulator.system->getSystemTime();
                simulator.simulateNumberOfSteps(intervalSteps, true);
                double endClock = simulator.system->getSystemTime();
                double elapsedTime = endClock - startClock;
                for (int i = 0; i < outputElectrodes.size(); ++i) {
                    int inCounts = 0;
                    int outCounts = 0;
                    double averageCurrent = 0.0;
                    for (int j = 0; j < numOfStates; ++j) {
                        outCounts += simulator.system->getNumberOfEvents(nAcceptors+outputElectrodes[i], j);
                        inCounts += simulator.system->getNumberOfEvents(j, nAcceptors+outputElectrodes[i]);            
                    } 
                    outputs[batch*outputElectrodes.size() + i] += static_cast<double>(inCounts-outCounts) / (elapsedTime*static_cast<double>(numOfIntervals));
                }
                simulator.system->resetEventCounts();
                intervalCounter++;
            }
        }
    }
    std::string fileName = saveFolderPath + "/batch" + std::to_string(batchID) + ".npz";
    cnpy::npz_save(fileName, "ID", &batchID, {1}, "w");
    cnpy::npz_save(fileName, "inputs", inputs.data(), shapeInputs, "a");
    cnpy::npz_save(fileName, "outputs", outputs.data(), shapeOutputs, "a");
}