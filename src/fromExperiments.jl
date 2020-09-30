using PyCall
using Plots; pyplot()
using Statistics
using LinearAlgebra

vtk = pyimport("vtk")
np = pyimport("numpy")

reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName("../data/h06.vtk")
reader.Update()
output = reader.GetOutput()
points = output.GetPoints()
r = [collect(points.GetPoint(i-1)) for i in 1:points.GetNumberOfPoints()]
length(r)

z = [r[i][3] for i=1:262144]

println(sum( abs.(z.-mean(z)))/length(z))
println(std(z))


z = reshape(z,512,512)
heatmap(z)

include("dfa.jl")
dfa2d(z)

function fbmArea(z)
    nx, ny = size(z)
    area = 0.0
    for i=1:nx-1
        for j=1:ny-1
            v1 = [0; 1; z[i,j+1] - z[i,j]]
            v2 = [1; 0; z[i+1,j] - z[i,j]]    
            v3 = [0; -1; z[i+1,j] - z[i+1,j+1]]
            v4 = [-1; 0; z[i,j+1] - z[i+1,j+1]]    
            area += (norm(v2 × v1)+norm(v4 × v3))/2
        end
    end
    
    Lp  = 0.0
    for i=1:nx-1 # eval Lp for each line
        
        for j=1:ny-1
            v1 = [1; 0; z[i+1,j] - z[i,j]]    
            Lp += norm(v1)
        end
    
    end
    
    return area/(nx*ny), (Lp/ny)/nx
end

fbmArea(z)

2^9


