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
	echo "Simulating device..."
	./main singleRun --equilibriumSteps=10000 --simulationSteps=10000000 --deviceName=test
}

build
run