import os

for i in range(1,51):
    os.system("julia Main.jl agis " + str(i))