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
	echo "Creating batch..."
	./main batchRun --batchSize=1000 --equilibriumSteps=10000 --simulationSteps=100000 --batchName=1
}

build
run