import os

os.system("julia ./MCF_min_congestion/Main.jl palmetto 30 ./MCF_min_congestion/results/")
os.system("julia ./MCF_min_congestion/Main.jl agis 30 ./MCF_min_congestion/results/")
os.system("julia ./MCF_min_congestion/Main.jl iowa 30 ./MCF_min_congestion/results/")
