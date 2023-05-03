import os

for i in range(1,21):
    os.system("julia makeInstances.jl /home/maxime/Téléchargements/MCF_energy/instances/agis/agis.graphml " + str(i) + " /home/maxime/Téléchargements/MCF_energy/instances/agis.nodetypes /home/maxime/Téléchargements/MCF_energy/instances/agis.nodecaps /home/maxime/Téléchargements/MCF_energy/instances/agis.edgecaps 3 [0.87,0.9,0.95] 1400 1800 20 40 150 200 0 2 8 5 15 /home/maxime/Téléchargements/MCF_energy/instances/agis/agis_bis" + str(i) +".jld")
    print(i)