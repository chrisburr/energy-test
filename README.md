# Optimised Energy Test

## Compiling

### macOS with gcc

This assumes gcc has been installed using homebrew.

```bash
g++-7 -o energy_test -std=c++1z -O3 -march=native -fopt-info -Wall -I. energy_test.cpp -fopenmp
```

### Manchester machines

```bash
source /cvmfs/sft.cern.ch/lcg/external/gcc/6.2.0/x86_64-slc6-gcc62-opt/setup.sh
g++ -o energy_test -std=c++1z -O3 -march=native -fopt-info -Wall -I. energy_test.cpp -fopenmp
```

## Input data

The data is assumed to be in space separated text files, see the example data in `data/csv/`. A script is included can create a suitable file (assuming `root_pandas` in installed):

```bash
./prepare_data.py  --input-fn=data/raw/sample0.root --output-fn-D0=data/csv/sample0-d0.csv --output-fn-D0bar=data/csv/sample0-d0bar.csv
```

## Running

Run the energy test for the input data only:

```bash
time ./energy_test data/csv/sample0-d0bar.csv data/csv/sample0-d0.csv
```

Run the energy test for the input data only, limiting the input files to 10000 events:

```bash
time ./energy_test --max-events=10000 data/csv/sample0-d0bar.csv data/csv/sample0-d0.csv
```

Run the energy test with permutations for 10000 D0 events and 11000 D0bar events:

```bash
time ./energy_test --max-events-1=10000 --max-events-2=11000 --n-permutations=100 data/csv/sample0-d0bar.csv data/csv/sample0-d0.csv
```
