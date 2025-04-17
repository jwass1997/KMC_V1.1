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