# Energy test with optimised permutations

[![DOI](https://zenodo.org/badge/113197119.svg)](https://zenodo.org/badge/latestdoi/113197119)

The [energy test](https://arxiv.org/abs/math/0309164) is a method for determining if two samples orginate from the same underlying distribution. This implimention can efficiently make use of multiple CPUs, with support for using scaled permutations as described in [TODO].

For an implimentation that can be used with NVIDIA GPUs, which may be faster for unscaled calculations with large samples (> 10⁷ points), see [Manet](https://manet.hepforge.org/).

## Compiling

### macOS with gcc

This assumes gcc has been installed using homebrew.

```bash
g++-7 -o energy_test -std=c++1z -O3 -march=native -Wall -I. energy_test.cpp -fopenmp
```

### Manchester machines

```bash
source /cvmfs/sft.cern.ch/lcg/external/gcc/6.2.0/x86_64-slc6-gcc62-opt/setup.sh
g++ -o energy_test -std=c++1z -O3 -march=native -Wall -I. energy_test.cpp -fopenmp
```

### lxplus

```bash
source /cvmfs/sft.cern.ch/lcg/external/gcc/6.2.0/x86_64-slc6-gcc62-opt/setup.sh
g++ -o energy_test -std=c++1z -O3 -Wall -I. energy_test.cpp -fopenmp
```

## Input data

The data is assumed to be in space separated text files and a script is included can create a suitable file (assuming `root_pandas` in installed):

```bash
./prepare_data.py  --input-fn=data/raw/sample0.root --output-fn-D0=data/csv/sample0-d0.csv --output-fn-D0bar=data/csv/sample0-d0bar.csv
```

## Running

Run the energy test for the input data only:

```bash
time ./energy_test sample1.csv sample2.csv
```

Run the energy test for the input data only, limiting the input files to 10000 events:

```bash
time ./energy_test --max-events=10000 sample1.csv sample2.csv
```

Run the energy test with permutations for 10000 points from sample 1 and 11000 points from sample 2:

```bash
time ./energy_test --max-events-1=10000 --max-events-2=11000 --n-permutations=100 sample1.csv sample2.csv
```
