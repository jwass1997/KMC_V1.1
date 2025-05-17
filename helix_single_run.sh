#!/bin/bash

ROOT_DIR="$(cd "$(dirname "$0")"/.. && pwd)"
SRC_DIR="$ROOT_DIR/KMC_V1.1"
LOCAL_DIR="$ROOT_DIR/local"

build() {
    echo "Compiling…"

    g++ -std=c++17 -Ofast \
        -I"$LOCAL_DIR/include"                 \
        "$SRC_DIR"/FEMmethods.cpp \
        "$SRC_DIR"/SystemGraph.cpp \
        "$SRC_DIR"/Simulator.cpp \
        "$SRC_DIR"/utils.cpp \
        "$SRC_DIR"/main.cpp \
        -I"${BOOST_INC_DIR}" \
        -L"$LOCAL_DIR/lib"                     \
        -L"${BOOST_LIB_DIR}" \
        -Wl,-rpath,"\$ORIGIN/../local/lib"     \
        -lmfem -lm -lcnpy -lz -lboost_program_options -fopenmp \
        -o "$SRC_DIR"/main || { 
            echo "Compilation FAILED" >&2
            exit 1
        }

    echo "Compiled to $SRC_DIR/main"
}

run() {
    echo "Running…"
    cd "$SRC_DIR"
    ./main singleRun \
         --equilibriumSteps=100000 \
         --simulationSteps=10000000 \
         --deviceName=test
}

build
run