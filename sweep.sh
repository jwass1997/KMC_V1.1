#!/usr/bin/env bash

build() {
	echo "Compiling..."
	g++ -std=c++17 -O3 \
	FEMmethods.cpp SystemGraph.cpp Simulator.cpp utils.cpp main.cpp \
	-o main \
	-lmfem -lm -lcnpy -lz -fopenmp -lboost_program_options \
	|| { echo "Compilation failed"; exit 1; }
}

run() {
	echo "Starting to sweep..."
	./main voltageSweep --numOfPoints=50 --equilibriumSteps=10000 --simulationSteps=100000 --fileName=test --min=-2.0 --max=2.0
}

build
run