Terminal commands and examples for running trace analysis 

Compiling the analyzer:
wsl g++ -O2 -std=c++17 -o trace_analyzer trace_analyzer.cpp -lz 

Running analyzer and exporting results: 
wsl ./trace_analyzer /mnt/c/Users/swe20/cbp7/trace_files/cbp-ng_training_traces/python3-pyperf-mako_2447_trace.gz > mako_results.txt 


Compiling the decoder: 
wsl g++ -O2 -std=c++17 -o trace_decode trace_decoder.cpp -lz

# default 100k
wsl ./trace_decode /mnt/c/Users/swe20/cbp7/trace_files/cbp-ng_training_traces/python3-pyperf-mako_2447_trace.gz mako_decoded.txt

# custom limit
wsl ./trace_decode /mnt/c/Users/swe20/cbp7/trace_files/cbp-ng_training_traces/python3-pyperf-mako_2447_trace.gz mako_decoded.txt 500000

Viewing exported results in VS code:
code mako_decoded.txt