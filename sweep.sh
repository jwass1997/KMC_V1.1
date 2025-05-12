#!/usr/bin/env bash

build() {
	echo "Compiling..."
	g++ -std=c++17 -Ofast \
	FEMmethods.cpp SystemGraph.cpp Simulator.cpp utils.cpp main.cpp \
	-o main \
	-lmfem -lm -lcnpy -lz -fopenmp -lboost_program_options \
	|| { echo "Compilation failed"; exit 1; }
}

run() {
	echo "Starting to sweep..."
	./main voltageSweep --numOfPoints=100 --equilibriumSteps=10000 --simulationSteps=100000 --fileName=warmUp1e4_sim1e6_sig0.0_2 --min=-1.5 --max=1.5
}

build
run