# Optimised Energy Test

## Compiling

### macOS with gcc

This assumes gcc has been installed using homebrew.

```bash
g++-7 -o energy_test -std=c++1z -Ofast -march=native -fopt-info -Wall -I. energy_test.cpp -fopenmp
```

### Manchester machines

```bash
source /cvmfs/sft.cern.ch/lcg/external/gcc/6.2.0/x86_64-slc6-gcc62-opt/setup.sh
g++ -o energy_test -std=c++1z -Ofast -march=native -fopt-info -Wall -I. energy_test.cpp -fopenmp
```

## Running

Prepare the data:

```bash
./prepare_data.py  --input-fn=data/raw/sample0.root --output-fn-D0=data/csv/sample0-d0.csv --output-fn-D0bar=data/csv/sample0-d0bar.csv
```

Run the energy test:

```bash
time ./energy_test --max-events=10000 --n-permutations=100 data/csv/sample0-d0bar.csv data/csv/sample0-d0.csv
```

