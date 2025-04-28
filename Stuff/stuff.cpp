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

    std::string fileName = saveFolderPath + "batch" + std::to_string(batchID) + ".npz";

    cnpy::npz_save(fileName, "ID", &batchID, {1}, "w");
    std::vector<double> inputs(batchSize*inputElectrodes.size(), 0.0);
    std::vector<double> outputs(batchSize*outputElectrodes.size(), 0.0);   
    std::vector<size_t> shapeInputs = {static_cast<size_t>(batchSize), inputElectrodes.size()}; 
    std::vector<size_t> shapeOutputs = {static_cast<size_t>(batchSize), outputElectrodes.size()}; 

    for (int batch = 0; batch < batchSize; ++batch) {
        Simulator simulator(defaultConfigs);
        int nAcceptors = simulator.system->getAcceptorNumber();
        int numOfStates = simulator.system->getNumOfStates();

        for (int inputElectrodeIndex = 0; inputElectrodeIndex < inputElectrodes.size(); ++inputElectrodeIndex) {
            double voltage = sampleFromUniformDistribution(minVoltage, maxVoltage);
            simulator.system->updateElectrodeVoltage(inputElectrodeIndex, voltage);
            inputs[inputElectrodeIndex] = voltage;
        }

        for (const auto& outputElectrodeIndex : outputElectrodes) {
            simulator.system->updateElectrodeVoltage(outputElectrodeIndex, 0.0);
        }

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
                outputs[i + batch*outputElectrodes.size()] += static_cast<double>(inCounts-outCounts) / (elapsedTime*static_cast<double>(numOfIntervals));
            }
            simulator.system->resetEventCounts();
            intervalCounter++;
        }
    }

    cnpy::npz_save(fileName, "inputs", inputs.data(), shapeInputs, "a");
    cnpy::npz_save(fileName, "outputs", outputs.data(), shapeOutputs, "a");
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

void voltageCurrentCharacteristic(
    double minVoltage, 
    double maxVoltage, 
    int numOfPoints, 
    int controlElectrodeIndex, 
    int scanElectrodeIndex, 
    int equilibriumSteps, 
    int simulationSteps, 
    int ID,
    const std::string& defaultConfig,
    const std::string& saveFolderPath) {

    if (saveFolderPath.empty()) {
        throw std::invalid_argument("No save folder specified !");
    }
    std::string dataFileName = saveFolderPath + "/IVcharacteristic" + std::to_string(ID) + ".npz";
    cnpy::npz_save(dataFileName, "ID", &ID, {1}, "w");

    double voltageStep = (maxVoltage - minVoltage) / (numOfPoints-1);

    std::vector<double> currentValues(numOfPoints, 0.0);
    std::vector<double> voltageValues(numOfPoints, 0.0);

    for (int i = 0; i < numOfPoints; ++i) {
        voltageValues[i] = minVoltage + voltageStep*i;
    }
    
    Simulator simulator(defaultConfig);

    for (int scanPoint = 0; scanPoint < numOfPoints; ++scanPoint) {
        std::cout << scanPoint << "\n";
        double currentValue = calculateCurrentAverage(
            simulator,
            controlElectrodeIndex,
            scanElectrodeIndex,
            voltageValues[scanPoint],
            equilibriumSteps,
            simulationSteps,
            100
        );
        currentValues[scanPoint] = currentValue;
    }
    
    std::vector<size_t> shape = {static_cast<size_t>(currentValues.size())};
    cnpy::npz_save(dataFileName, "currentValues", currentValues.data(), shape, "a");
}

double calculateCurrentAverage(
    Simulator& simulator, 
    int controlElectrodeIndex,
    int scanElectrodeIndex, 
    double voltage, 
    int equilibriumSteps, 
    int simulationSteps, 
    int numberOfIntervals) {

    /**
     * 
     * Single current average for a specific voltage
     * 
     */
    simulator.system->updateElectrodeVoltage(controlElectrodeIndex, voltage);
    simulator.system->updateElectrodeVoltage(scanElectrodeIndex, 0.0);
    simulator.simulateNumberOfSteps(equilibriumSteps, false);
    
    int nAcceptors = simulator.system->getAcceptorNumber();
    int numOfStates = simulator.system->getNumOfStates();

    double averageCurrent = 0.0;
    int intervalSteps = simulationSteps / numberOfIntervals;
    int intervalCounter = 0;
    while (intervalCounter < numberOfIntervals) {
        double startClock = simulator.system->getSystemTime();
        simulator.simulateNumberOfSteps(intervalSteps, true);
        double endClock = simulator.system->getSystemTime();

        int inCounts = 0;
        int outCounts = 0;
        for (int i = 0; i < numOfStates; ++i) {
            outCounts += simulator.system->getNumberOfEvents(nAcceptors+scanElectrodeIndex, i);
            inCounts += simulator.system->getNumberOfEvents(i, nAcceptors+scanElectrodeIndex);            
        }
        double elapsedTime = endClock - startClock;
        averageCurrent += static_cast<double>(inCounts-outCounts) / (elapsedTime*static_cast<double>(numberOfIntervals));

        simulator.system->resetEventCounts(); 
        intervalCounter++;
    }

    return averageCurrent;
}

/**
 * 
 * Old main function
 * 
 */
std::string saveFolder = "currentData";
    std::string defaultDeviceConfigs = "default_configs";
    createDirectoryFromStringPath("", saveFolder);
    int controlElectrodeIndex = 1;
    int scanElectrodeIndex = 0;

    std::string deviceID1 = "1";
    std::vector<int> inputElectrodes = {1,2,3,4,5,6,7};
    std::vector<int> outpuElectrodes = {0};
    int maxThreads = omp_get_max_threads();
    //std::cout << maxThreads << "\n";
    //std::cout<<__cplusplus<<"\n";

    int numOfPoints = 100;
    double maxVoltage = -1.;
    double minVoltage = 1.;
    double range = maxVoltage - minVoltage;
    double voltageStep = range / (numOfPoints-1);
    std::vector<double> currents(numOfPoints, 0.0);
    #pragma omp parallel for
    for (int voltagePoint = 0; voltagePoint < numOfPoints; ++voltagePoint) {
        std::cout<<voltagePoint<<"\n";
        std::vector<double> voltageSetting = {0.,0.,0.,0.,0.,0.,0.,0.};
        voltageSetting[1] = minVoltage + voltageStep*voltagePoint;
        double outputCurrent = 0.0;
        outputCurrent = IVPoint(voltageSetting, 5, 0, 1e4, 1e7, 100, voltagePoint, defaultDeviceConfigs, saveFolder);
        currents[voltagePoint] = outputCurrent;
    }
    std::string fileName = "currents.npz";
    cnpy::npz_save(fileName, "currents", currents.data(), {static_cast<size_t>(numOfPoints)}, "w");