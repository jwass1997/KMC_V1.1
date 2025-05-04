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
            stateEnergies[nAcceptors + electrodeIndices[i]] = newVoltages[i]*e / kbT;
        }
    
        resetPotential();
    
        finiteElementSolver->updateMultipleElectrodes(electrodeIndices, newVoltages);
        finiteElementSolver->solvePoissonEquation();
    
        updatePotential();
    }

/* else if (firstCommand == "IVpoint") {
        boost::program_options::options_description IVRunOptions("IV point calculation");
        IVRunOptions.add_options()
            ("voltage", boost::program_options::value<double>()->required())
            ("numOfDevices", boost::program_options::value<int>()->default_value(30))
            ("scanElectrodeIndex", boost::program_options::value<int>()->default_value(0))
            ("equilibriumSteps", boost::program_options::value<int>()->default_value(1e4))
            ("simulationSteps", boost::program_options::value<int>()->required())
            ("numOfIntervals", boost::program_options::value<int>()->default_value(100))
            ("ID", boost::program_options::value<int>()->required())
            ("saveFolderPath", boost::program_options::value<std::string>()->required())
        ;
            
        boost::program_options::variables_map IVRunVM;
        boost::program_options::store(
            boost::program_options::command_line_parser(
                remainingCommand).options(IVRunOptions).run(),
                IVRunVM);
        boost::program_options::notify(IVRunVM);

        std::vector<double> voltageSetting(8, 0.0);
        voltageSetting[(IVRunVM["scanElectrodeIndex"].as<int>() + 1) % 8] = IVRunVM["voltage"].as<double>();

        double current = IVPoint(
            voltageSetting,
            IVRunVM["numOfDevices"].as<int>(),
            IVRunVM["scanElectrodeIndex"].as<int>(),
            IVRunVM["equilibriumSteps"].as<int>(),
            IVRunVM["simulationSteps"].as<int>(),
            IVRunVM["numOfIntervals"].as<int>(),
            IVRunVM["ID"].as<int>(),
            "default_configs"
        );

        std::string fileName = IVRunVM["saveFolderPath"].as<std::string>() + "/point_" + std::to_string(IVRunVM["ID"].as<int>()) + ".npz";
        cnpy::npz_save(fileName, "input", &IVRunVM["voltage"].as<double>(), {1}, "w");
        cnpy::npz_save(fileName, "output", &current, {1}, "a");

        return 1;
    } */

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
    int numPoints = 50;
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
        -0.7,
        1.1,
        -0.8,
        0.7,
        -1.2,
        1.0
    };
    sim.system->updateVoltages(voltageSetting);
    /* FiniteElementeBase* finElem = sim.system->finiteElementSolver;

    mfem::GridFunction potential = *(sim.system->finiteElementSolver->solutionVector);
    auto mesh = finElem->mesh;
    auto sol  = finElem->solutionVector;
    int maxNumberOfElements = 1e5;
    std::ofstream out("potential.csv");
    out << "x,y,potential\n";
    for (int i = 0; i < mesh->GetNV(); i++)
    {
        auto vert = mesh->GetVertex(i);             // a double[2]
        double phi = (*sol)[i];                     // solution at that DOF
        out << vert[0] << "," << vert[1] << "," << phi << "\n";
    }
    out.close(); */
    
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

    /* for (int i = 0; i < nS; ++i) {
        std::cout << sim.system->stateEnergies[i] << "\n";
    } */

    /* for (int i = 0; i < nAcceptors; ++i) {
        std::cout << sim.system->acceptorCoordinates[i*2] << " " << sim.system->acceptorCoordinates[i*2+1] << "\n";
    } */

    std::string defaultConfig = "default_configs";

    Simulator sim(defaultConfig);

    int nA = sim.system->nAcceptors;
    int nE = sim.system->nElectrodes;
    int numOfStates = 200 + nE;
    int numPoints = 200;
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
        -0.7,
        1.1,
        -0.8,
        0.7,
        -1.2,
        1.0
    };
    sim.system->updateVoltages(voltageSetting);

    int outputElectrodeIdx = 0;
    for (int i = 0; i < numPoints; ++i) {
        voltageSetting[1] = minVoltage + i*vStep;
        sim.system->updateVoltages(voltageSetting);
        current = 0.0;
        sim.simulateNumberOfSteps(1e4, false);
        sim.system->resetEventCounts();
        double startTime = sim.system->getSystemTime();
        sim.simulateNumberOfSteps(1e4, true);
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
    file.close();

    double myGetPotential(const double& x, const double& y, mfem::GridFunction solutionVector, int maxNumberOfElements) {
        double radius = std::sqrt(2.0*M_PI*150.0*150.0 / 200.0);
        int layers = std::sqrt(maxNumberOfElements / 6.0);
        double deltaR = radius / layers;
        int layer = std::sqrt(x * x + y * y) / deltaR + 0.5;
        double phi = std::atan2(y, x) / (2 * PI);
        // std::cout<<"return index: "<<int(3*(layers+1)*layers-3*layer*(layer+1) +
        // (phi < 0 ? phi + 1 : phi) * 6 * layer + 0.5)<<std::endl;
        return solutionVector[int(3 * (layers + 1) * layers - 3 * layer * (layer + 1) + (phi < 0 ? phi + 1 : phi) * 6 * layer + 0.5)];
    } 