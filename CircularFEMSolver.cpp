#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#include "CircularFEMSolver.h"

CircularFEMSolver::CircularFEMSolver(
    double radius, 
    int numOfCircleVertices, 
    int numOfConcentricCircles, 
    std::vector<Electrode* >& electrodes
)
    : radius(radius)
    , numOfCircleVertices(numOfCircleVertices)
    , numOfConcentricCircles(numOfConcentricCircles)
    , mesh(nullptr)
    , solution(nullptr)
    , fec(nullptr)
    , fespace(nullptr)
    , a(nullptr)
    , b(nullptr)
{
    if (numOfCircleVertices < 3 || numOfConcentricCircles < 1) {
        throw std::invalid_argument("Circular mesh needs at least 3 vertices per circle and at least one circle");
    }

    electrodePars = electrodes;
    electrodesTdofIndices.resize(electrodePars.size());

    createCircularMesh(numOfCircleVertices, numOfConcentricCircles);
    if (order > 0) {
        fec = new H1_FECollection(order, dim);
    }
    fespace = new FiniteElementSpace(mesh, fec);
    solution = new GridFunction(fespace);
    *solution = 0;

    setElectrodes();
}

CircularFEMSolver::~CircularFEMSolver() {
    delete mesh;
    delete solution;
    delete fec;
    delete fespace;
}

void CircularFEMSolver::createCircularMesh(int numOfCircleVertices, int numOfConcentricCircles) {

    int numOfVertices = 1 + numOfConcentricCircles*numOfCircleVertices;    
    int numOfElements = numOfCircleVertices + (numOfConcentricCircles-1)*2*numOfCircleVertices;
    int numOfBdrElements = numOfCircleVertices;
    mesh = new Mesh(dim, numOfVertices, numOfElements, numOfBdrElements);    

    circleVertices.SetSize(numOfVertices);
    circleVertices[0](0) = 0.0;
    circleVertices[0](1) = 0.0;
    mesh->AddVertex(circleVertices[0](0), circleVertices[0](1));

    for (int i = 1; i < numOfConcentricCircles+1; ++i) {
        double circleRadius = radius*static_cast<double>(i) / static_cast<double>(numOfConcentricCircles);
        for (int j = 1; j < numOfCircleVertices+1; ++j) {
            double angle = 2.0*M_PI*static_cast<double>(j) / static_cast<double>(numOfCircleVertices);
            
            double xCoord = circleRadius*std::cos(angle);
            double yCoord = circleRadius*std::sin(angle);

            circleVertices[1 + (i-1)*numOfCircleVertices + (j-1)](0) = xCoord;
            circleVertices[1 + (i-1)*numOfCircleVertices + (j-1)](1) = yCoord;
            mesh->AddVertex(circleVertices[1 + (i-1)*numOfCircleVertices + (j-1)](0), circleVertices[1 + (i-1)*numOfCircleVertices + (j-1)](1));
        }        
    }

    int innerTriangle[3];
    for (int i = 1; i < numOfCircleVertices+1; ++i) {
        int vertexIndex1 = i;
        int vertexIndex2 = i % numOfCircleVertices + 1;
        innerTriangle[0] = 0;
        innerTriangle[1] = vertexIndex1;
        innerTriangle[2] = vertexIndex2;
        mesh->AddTriangle(innerTriangle);
    }

    int outerTriangle1[3];
    int outerTriangle2[3];
    for (int i = 1; i < numOfConcentricCircles; ++i) {
        for (int j = 1; j < numOfCircleVertices+1; ++j) {
            outerTriangle1[0] = (i-1)*numOfCircleVertices + j;
            outerTriangle1[1] = i*numOfCircleVertices + j;
            outerTriangle1[2] = i*numOfCircleVertices + j % numOfCircleVertices + 1;

            outerTriangle2[0] = (i-1)*numOfCircleVertices + j;
            outerTriangle2[1] = (i-1)*numOfCircleVertices + j % numOfCircleVertices + 1;
            outerTriangle2[2] = i*numOfCircleVertices + j % numOfCircleVertices + 1;

            mesh->AddTriangle(outerTriangle1);
            mesh->AddTriangle(outerTriangle2);
        }
    }

    int bdrElement[2];
    for (int i = 0; i < numOfCircleVertices; ++i) {
        bdrElement[0] = 1 + (numOfConcentricCircles-1)*numOfCircleVertices + i;
        bdrElement[1] = 1 + (numOfConcentricCircles-1)*numOfCircleVertices + (i+1) % numOfCircleVertices;
        mesh->AddBdrSegment(bdrElement, 1);
    }
    mfem::out.Disable();
    mesh->FinalizeTriMesh(0, 0, true);
}

void CircularFEMSolver::solvePoissonEquation() {

    a = new BilinearForm(fespace);
    ConstantCoefficient one(1.0);
    a->AddDomainIntegrator(new DiffusionIntegrator(one));
    a->Assemble();

    b = new LinearForm(fespace);
    ConstantCoefficient zero(0.0);
    b->AddDomainIntegrator(new DomainLFIntegrator(zero));
    b->Assemble();

    a->FormLinearSystem(ess_tdof_list, *solution, *b, A, X, B);

    SparseMatrix &A_sp = *(dynamic_cast<SparseMatrix*>(A.Ptr()));
    M = GSSmoother(A_sp);
    PCG(A_sp, M, B, X, 0, 1000, 1e-12, 0.0); 
    a->RecoverFEMSolution(X, *b, *solution);

    delete a;
    delete b;
}

double CircularFEMSolver::calculatePotential(double x, double y) {

    double r = std::sqrt(x*x + y*y);
    double phi = std::atan2(y, x);

    if (phi < 0.0) {
        phi += 2.0*M_PI;
    }

    double dr = radius / static_cast<double>(numOfConcentricCircles);
    double dphi = 2.0*M_PI / static_cast<double>(numOfCircleVertices);

    int circleIndex = std::min(numOfConcentricCircles-1, static_cast<int>(r/dr));
    int vertexIndex = int(std::round(phi/dphi)) % numOfCircleVertices;

    return (*solution)[1 + circleIndex*numOfCircleVertices + vertexIndex];
}

void CircularFEMSolver::updateElectrode(double angularPosition, int electrodeIndex, double newVoltage) {
    
    //electrodePars[electrodeIndex]->voltage = newVoltage;

    std::vector<int> indices = getElectrodeVertices(angularPosition);

    for (const auto& index : electrodesTdofIndices[electrodeIndex]) {
        (*solution)[ess_tdof_list[index]] = newVoltage;
    }
}

void CircularFEMSolver::updateMultipleElectrodes(std::vector<int> electrodeIndices, std::vector<double> newVoltages) {

    for (int i = 0; i < electrodeIndices.size(); ++i) {
        for (const auto& index : electrodesTdofIndices[electrodeIndices[i]]) {
            (*solution)[ess_tdof_list[index]] = newVoltages[i];
        }
    }
}

std::vector<int> CircularFEMSolver::getElectrodeVertices(double angularPosition) {

    if (angularPosition < 0.0 || angularPosition > 360.0) {
        throw std::invalid_argument("Angular position not normalized");
    }

    double deltaPhi = 1.0 / 20.0;
    angularPosition = angularPosition / 360.0;

    bool includingZero = false;

    double start = std::floor((angularPosition - deltaPhi)*static_cast<double>(numOfCircleVertices));
    double end = std::ceil((angularPosition + deltaPhi)*static_cast<double>(numOfCircleVertices));

    if (start < 0) {
        includingZero = true;
        start += numOfCircleVertices;
    }

    start += (numOfConcentricCircles-1)*numOfCircleVertices + 1;
    end += (numOfConcentricCircles-1)*numOfCircleVertices + 1;

    std::vector<int> electrodeVertices;
    int nBdr = mesh->GetNBE();

    for (int i = 0; i < nBdr; ++i) {
        Element* bdrElement = mesh->GetBdrElement(i);
        int* bdrVertices = bdrElement->GetVertices();        
        int numOfBdrvertices = bdrElement->GetNVertices();
        for (int j = 0; j < 2; ++j) {
            if (!includingZero) {
                if (bdrVertices[j] >= start && bdrVertices[j] <= end) {
                    electrodeVertices.push_back(bdrVertices[j]);
                }
                else {
                    continue;
                }
            }
            if (includingZero) {
                if (bdrVertices[j] >= start || bdrVertices[j] <= end) {
                    electrodeVertices.push_back(bdrVertices[j]);
                }
                else {
                    continue;
                }
            }
        }
    }

    std::sort(electrodeVertices.begin(), electrodeVertices.end());
    auto it = std::unique(electrodeVertices.begin(), electrodeVertices.end());
    electrodeVertices.erase(it, electrodeVertices.end());

    return electrodeVertices;
}

void CircularFEMSolver::setElectrodes() {
    int tdofIndex = 0;
    int electrodeIndex = 0;
    for (const auto& electrode : electrodePars) {
        double angularPosition = electrode->angularPosition;
        double voltage = electrode->voltage;

        std::vector<int> electrodeVertices = getElectrodeVertices(angularPosition);
        for (const auto& index : electrodeVertices) {
            (*solution)[index] = voltage;
            ess_tdof_list.Append(index);
            electrodesTdofIndices[electrodeIndex].push_back(tdofIndex);
            tdofIndex++;
        }
        electrodeIndex++;
    }
}

void exportPotentialToCSV(CircularFEMSolver &solver, const std::string &filename, int numPoints, double R)
{
    std::ofstream ofs(filename);
    if (!ofs.is_open())
    {
        std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
        return;
    }

    // Write CSV header
    ofs << "x,y,potential\n";

    // Create a uniform grid over the square [-R, R] x [-R, R]
    for (int i = 0; i < numPoints; ++i)
    {
        // x coordinate ranging from -R to R
        double x = -R + (2 * R) * i / (numPoints - 1);
        for (int j = 0; j < numPoints; ++j)
        {
            // y coordinate ranging from -R to R
            double y = -R + (2 * R) * j / (numPoints - 1);
            // Only consider points within the circle
            if (x * x + y * y <= R * R)
            {
                double potential = solver.calculatePotential(x, y);
                ofs << x << "," << y << "," << potential << "\n";
            }
        }
    }
    ofs.close();
}