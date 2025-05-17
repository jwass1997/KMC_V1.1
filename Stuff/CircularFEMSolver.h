#pragma once

#include <vector>
#include <memory>

#include "mfem.hpp"
#include "utils.h"

using namespace mfem;

class CircularFEMSolver{
public:
    CircularFEMSolver();

    CircularFEMSolver(
        double radius, 
        int numOfCircleVertices,
        int numOfConcentricCircles, 
        std::vector<Electrode*> &electrodes
    );

    ~CircularFEMSolver();

    const int dim = 2;
    const int sdim = 2;
    const int order = 1;

    double radius = 1.0;
    int numOfCircleVertices;
    int numOfConcentricCircles;

    std::vector<Electrode* > electrodePars;

    Array<Vertex> circleVertices;

    Mesh* mesh;
    GridFunction* solution;
    FiniteElementCollection* fec;
    FiniteElementSpace* fespace;
    BilinearForm* a;
    LinearForm* b;
    OperatorPtr A;
    Vector B, X;
    GSSmoother M;
    Array<int> ess_tdof_list;

    std::vector<std::vector<int>> electrodesTdofIndices;

    void createCircularMesh(int numOfCircleVertices, int numOfConcentricCircles);

    void solvePoissonEquation();

    double calculatePotential(double x, double y);

    void updateElectrode(
        double angularPosition, 
        int electrodeIndex, 
        double newVoltage);

    void updateMultipleElectrodes(std::vector<int> electrodeIndices, std::vector<double> newVoltages);

    std::vector<int> getElectrodeVertices(double position);
    
    void setElectrodes();
};