# Dynamic-SAM-Research

Research repository for Dynamic Suffix Automaton / Suffix Tree algorithms and comparisons.

## Contents
- `prev_algo/` — previous algorithm implementations used for comparison:
  - `bruteforce+hashing.cpp`
  - `Dynamic Suffix Trees version 2.cpp`
  - `naive_bruteforce.cpp`
  - `Suffix Array + LCP + RMQ.cpp`
  - `Treap-Based APSPM.cpp`
- `genarator and dataset/` — dataset generation and sample data
- `info and result/` — info and result text files
- `pictureee/` — plots and CSV used for comparisons

## Build / Run (example)
You can compile individual C++ files with `g++` (MinGW on Windows or any modern g++):

```
g++ -O2 -std=c++17 "prev_algo\\bruteforce+hashing.cpp" -o bruteforce_hashing.exe
.\\bruteforce_hashing.exe
```

Replace the file path/name for other implementations in `prev_algo/` as needed.

## Notes
- This repository currently has no license file.
- If you want, add a `.gitignore` to exclude build artifacts and editor files.

## Contact
For questions about the code or datasets, open an issue or contact the maintainer.
