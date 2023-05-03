include("makeInstances.jl")

for i in 1:30
    main(split("/home/maxime/Téléchargements/MCF_energy/instances/agis/agis.graphml " * string(i) * " /home/maxime/Téléchargements/MCF_energy/instances/agis/agis.nodetypes /home/maxime/Téléchargements/MCF_energy/instances/agis/agis.nodecaps /home/maxime/Téléchargements/MCF_energy/instances/agis/agis.edgecaps 3 [0.87,0.9,0.95] 1400 1800 20 40 150 200 0 2 8 5 15 /home/maxime/Téléchargements/MCF_energy/instances/agis/agis" * string(i) *".jld"))
end

for i in 1:30
    main(split("/home/maxime/Téléchargements/MCF_energy/instances/iowa/iowa.graphml " * string(i) * " /home/maxime/Téléchargements/MCF_energy/instances/iowa/iowa.nodetypes /home/maxime/Téléchargements/MCF_energy/instances/iowa/iowa.nodecaps /home/maxime/Téléchargements/MCF_energy/instances/iowa/iowa.edgecaps 3 [0.87,0.9,0.95] 1400 1800 20 40 150 200 0 2 8 5 15 /home/maxime/Téléchargements/MCF_energy/instances/iowa/iowa" * string(i) *".jld"))
end

for i in 1:30
    main(split("/home/maxime/Téléchargements/MCF_energy/instances/palmetto/palmetto.graphml " * string(i) * " /home/maxime/Téléchargements/MCF_energy/instances/palmetto/palmetto.nodetypes /home/maxime/Téléchargements/MCF_energy/instances/palmetto/palmetto.nodecaps /home/maxime/Téléchargements/MCF_energy/instances/palmetto/palmetto.edgecaps 3 [0.87,0.9,0.95] 1300 2800 15 30 150 200 0 2 8 5 15 /home/maxime/Téléchargements/MCF_energy/instances/palmetto/palmetto" * string(i) *".jld"))
end