# BMF_LNS

**Boolean Matrix Factorisation with GI, LS, WLS, LNS**

## Repository

```
git clone https://github.com/Smartfreak14/BMF_LNS.git
cd BMF_LNS
```

## Description

A C++ solver for Boolean Matrix Factorization (BMF) using MaxSAT-based local search methods:
- **GI** (Greedy Initialize) -- greedy construction of initial factorization A, B
- **LS** (Local Search) -- 1-flip hill climbing on A and B
- **WLS** (Weighted Local Search) -- weighted variant with dynamic penalty
- **LNS** (Large Neighborhood Search) -- destroy-and-repair via partial MaxSAT

Given a binary matrix M, finds A (m x k) and B (k x n) such that the Boolean product A . B approximates M with minimal reconstruction errors.

## Build

```bash
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4
```

## Usage

```bash
# Single run with options (opt1 to opt6) move the file path to execute on all data
./rbac_main --opt1 ../data/Real/zoo.csv 25 30        # GREEDY + LS
./rbac_main --opt2 ../data/Real/zoo.csv 25 30        # GREEDY + Fix_A/Fix_B
./rbac_main --opt3 ../data/Real/zoo.csv 25 30        # GREEDY + LNS
./rbac_main --opt4 ../data/Real/zoo.csv 25 30        # GREEDY + LS + LNS
./rbac_main --opt5 ../data/Real/zoo.csv 25 30        # GREEDY + Fix_A/Fix_B + LNS
./rbac_main --opt6 ../data/Real/zoo.csv 25 30        # GREEDY + LS + Fix_A/Fix_B + LNS

# for helper option
./rbac_main --help

# With timeout (seconds)
./rbac_main --opt2 ../data/Real/apj.csv 454 60
```

Arguments: `<csv_file> <k> <timeout_seconds>`

## Project Structure

```
src/
  main.cpp              # Main program with opt1-opt6
  BMF.cpp / BMF.hpp     # MaxSAT-based BMF solver (GI, LNS)
  BMFLocalSearch.cpp/hpp # Local search (LS, WLS)
  SATSolver.cpp/hpp     # SAT/MaxSAT interface (EvalMaxSAT)
  Matrix.hpp            # Matrix data structure
  CSVMatrixLoader.hpp   # CSV reader
  DataUtils.hpp         # Configuration & results I/O
  verif.py              # Python verification script
data/
  Real/                 # Real-world benchmark datasets (31 CSV)
  Synthetic/            # Synthetic datasets (12 CSV)
third_party/
  EvalMaxSAT/           # MaxSAT solver
  minisat220/           # SAT solver backend
  MaLib/                # Utility library
```

## Dependencies

- C++17 compiler (GCC, Clang)
- CMake >= 3.10
- EvalMaxSAT (included in third_party/)

## License

See individual licenses in `third_party/` for EvalMaxSAT, CaDiCaL, Glucose, and MiniSat.