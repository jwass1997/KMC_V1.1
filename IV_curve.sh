#!/usr/bin/env bash

mkdir -p IVCurve

MAXJOBS=10
count=0

build() {
	echo "Compiling..."
	g++ -std=c++17 -Ofast \
	CircularFEMSolver.cpp SystemGraph.cpp Simulator.cpp utils.cpp main.cpp \
	-o main \
	-lmfem -lm -lcnpy -lz -fopenmp -lboost_program_options \
	|| { echo "Compilation failed"; exit 1; }
}

build

echo "Running..."

for x in $(< voltages.txt); do
    echo "$count"
    ./main IVpoint \
    --voltage="$x" \
    --numOfDevices=30 \
    --scanElectrodeIndex=0 \
    --equilibriumSteps=1000 \
    --simulationSteps=1000 \
    --numOfIntervals=100 \
    --ID="$count" \
    --saveFolderPath=IVCurve &

    while (( $(jobs -rp | wc -l) >= MAXJOBS )); do
        wait -n
    done

    ((count++))

done

wait

echo "All $((count-1)) jobs complete"