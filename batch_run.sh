#!/usr/bin/env bash

build() {
	echo "Compiling..."
	g++ -std=c++17 -Ofast \
	CircularFEMSolver.cpp SystemGraph.cpp Simulator.cpp utils.cpp main.cpp \
	-o main \
	-lmfem -lm -lcnpy -lz -fopenmp -lboost_program_options \
	|| { echo "Compilation failed"; exit 1; }
}

run() {
	echo "Creating batch..."
	./main batchRun --batchSize=20 --equilibriumSteps=100000 --simulationSteps=1000000 --batchName=test
}

build
run