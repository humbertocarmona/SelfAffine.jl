function showarray(arr)
    println(summary(arr))
    Base.print_matrix(IOContext(stdout, :limit => true), arr)
    println()
end

function strided(carr::Array{Float64,1}, scale::Int64, offset=-1)
    offset = offset==-1 ? scale : offset
    stridedarr = [carr[n:n+scale-1] for n=1:offset:(length(carr)-scale+1)]
end

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
    return area/(nx*ny)
end
