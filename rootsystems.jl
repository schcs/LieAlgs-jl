using LinearAlgebra 

function simple_root_system( type::String, rank::Int64 )
    
    if !( type in ["A","D"] ) 
        throw( "Error: type $type is not yet implemented" )
    end

    simple_system = []

    if type == "A" 

        for i in 1:rank
            vec = zeros( Int64, rank+1 )
            vec[i], vec[i+1] = 1, -1
            push!( simple_system, vec )
        end

    elseif type == "D"

        for i in 1:rank-1
            vec = zeros( Int64, rank )
            vec[i], vec[i+1] = 1, -1
            push!( simple_system, vec )
        end
        vec = zeros( Int64, rank )
        vec[rank-1], vec[rank] = 1, 1
        push!( simple_system, vec )

    end 

    return simple_system
end

function cartan_matrix( simple_system::Vector )::Matrix
    
    l = length( simple_system )
    mat = zeros( Int64, l, l )
    for i in 1:l
        for j in 1:l
            mat[i,j] = 2*dot(simple_system[i],simple_system[j])/
                        dot(simple_system[j],simple_system[j])
        end
    end

    return mat
end

function cartan_matrix( m::Int64, n::Int64 )
    mat = zeros( Int64, m+n-1, m+n-1 )
    for i in 1:m
        mat[i,i] = 1
    end

    for i in m+1:m+n-1
        mat[i,i] = -1
    end

    return mat 
end


function complete_root_system( simple_system::Vector )

    l = length( simple_system )
    C = cartan_matrix( simple_system )

    phi = [ zeros( Int64, l ) for i in 1:l ]
    for i in 1:l 
        phi[i][i] = 1 
    end 

    start_wn = [1,l+1]
    n = 1
    enlarged = true
    while enlarged
        push!( start_wn, start_wn[end] ) 
        enlarged = false
        for i in start_wn[n]:(start_wn[n+1]-1)
            for j in 1:l
                gamma = Vector(copy( phi[i] ))

                r = gamma[j]
                gamma[j] = 0
                while !( gamma in phi )
                    gamma[j] += 1
                    r -= 1
                end
                q = r - sum( phi[i][u]*C[u,j] for u in 1:l )
                new_r = phi[i] + phi[j] 
                if q > 0 && ! (new_r in phi)
                    push!( phi, phi[i]+phi[j] )
                    enlarged = true
                    start_wn[end] += 1                     
                end 
            end
        end
        n += 1
    end
    return cat( phi, .-(phi), dims = 1 )
end

