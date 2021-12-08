# the type for structure constants tables
sc_type = Matrix{Vector{Tuple{Int64,Int64}}}


# the type for holding information about Lie algebras
struct lie_algebra{T}                     # T is the coefficient domain (Int64, etc)
    type::String                          # type is "A", "B", etc
    rank::Int64                           # Lie rank
    cartan_matrix::Matrix{Int64}          # Cartan matrix
    simple_roots::Vector{Vector{Int64}}   # list of simple roots
    root_system::Vector{Vector{Int64}}    # list of roots
    structure_constants::sc_type          # the structure constants table 
end

# the type for a Lie algebra element 
mutable struct lie_algebra_element{T}
    parent::lie_algebra                   # the parent Lie algebra
    content::Vector{T}                    # the element
end

# the function for constructing Lie algebras associated with root systems
# and for calculating with their elements
function lie_algebra( type::String, rank::Int64; dom = Int64 )::lie_algebra
    sr = simple_root_system( type, rank )
    cm = cartan_matrix( sr )
    root_system = sort( complete_root_system( sr ))
    
    dim = rank+length(root_system)
    table = fill( [], dim, dim )
    for i in 1:dim
        for j in i+1:dim
            m = mult( sr, root_system,  cm, i, j )
            table[i,j] = m
            table[j,i] = [ (-x[1],x[2]) for x in m ]
        end 
    end 
           
    return lie_algebra{dom}( type, rank, cm, sr, root_system, table )
end

function dimension( l::lie_algebra )

    return size( l.structure_constants )[1]
end 

function simple_roots( l::lie_algebra )
    return l.simple_roots
end 

function root_system( l::lie_algebra )
    return l.root_system
end 


function cartan_matrix( l::lie_algebra )
    return l.cartan_matrix
end 

function rank( l::lie_algebra )
    return l.rank 
end 

function lie_mon_to_string( alg, x )

    return "$(x[1])*$( x[2]<=rank(alg) ? "h($(-x[2]))" : "x($(alg.root_system[x[2]-rank(alg)]))")" 
end

function domain( l::lie_algebra )
    return typeof( l ).parameters[1]
end

function Base.display( x::lie_algebra_element )
    
    iszero = true
    for i in 1:dimension( x.parent )
        if x.content[i] != 0 
            println( lie_mon_to_string( x.parent, ( x.content[i], i )))
            iszero = false 
        end 
    end 

    if iszero 
        println( "0" )
    end 
end

function Base.display( l::lie_algebra )
    println( "Lie algebra of type $(l.type) and rank $(l.rank)" )
end 

function zero( l::lie_algebra )
    T = domain( l )
    return lie_algebra_element{T}( l, zeros( T, dimension( l )))
end

# The implementation of the N function to construct structure constants. Only in the 
# simpy laced case. See page 191 in de Graaf. 

function N( C, x, y )
    l = length( x )
    eps(i,j) = i == j ? -1 : (C[i,j]*C[j,i] != 0 && i<j ? -1 : 1 )
    return prod( eps(i,j)^(x[i]*y[j]) for i in 1:l for j in 1:l )
end   

# adds x and y and puts the result into x
# assumes x in ordered form

function add!( x, y )

    x.content += y.content
end 

function Base.:+(  x::lie_algebra_element, y::lie_algebra_element )

    return lie_algebra_element{domain( x.parent )}( x.parent, x.content+y.content )
end

function Base.:*(  x::Number, y::lie_algebra_element )

    if ! isa( x, domain( y.parent ))
        throw( "coefficient not of the right type ")
    end 
    return lie_algebra_element{domain( y.parent )}( y.parent, x*y.content )
end

function Base.:-(  y::lie_algebra_element )

    return lie_algebra_element{domain( y.parent )}( y.parent, -y.content )
end

function Base.:-(  x::lie_algebra_element, y::lie_algebra_element )

    return x+(-1)*y
end
    

function mult( simple_roots, roots, cart_mat, x, y )

    nr_simple_roots = length( simple_roots )
    nr_roots = length( roots )

    if x <= nr_simple_roots && y <= nr_simple_roots 
        return []
    elseif x <= nr_simple_roots
        ry = roots[y-nr_simple_roots]
        root = sum( ry[i]*simple_roots[i] for i in 1:nr_simple_roots )
        coeff = dot( root, simple_roots[x] ) 
        return coeff == 0 ? [] : [( coeff, y )]
    elseif y <= nr_simple_roots
        rx = roots[y-nr_simple_roots]
        root = sum( rx[i]*simple_roots[i] for i in 1:nr_simple_roots )
        coeff = dot( root, simple_roots[y] ) 
        return coeff == 0 ? [] : [( -coeff, y )]
    elseif roots[y-nr_simple_roots] == -roots[x-nr_simple_roots] 
        res = [ (-roots[x-nr_simple_roots][i],i) for i in 1:nr_simple_roots ]
        return [ x for x in res if x[1] != 0 ]
    else
        z = roots[x-nr_simple_roots]+roots[y-nr_simple_roots]
        pos = searchsortedfirst( roots, z, lt = (u,v) -> u<v )
        if pos <= nr_roots && roots[pos] == z
            res = [( N( cart_mat, roots[x-nr_simple_roots], roots[y-nr_simple_roots] ), 
                        pos+nr_simple_roots )]
            return [ x for x in res if x[1] != 0 ] 
        else 
            return []
        end 
    end
end 

#=
function Base.:*( x::lie_algebra_element, y::lie_algebra_element )

    if x.parent != y.parent 
        throw( "error: not elements of the same algebra" )
    end 
    res = []

    for i in x.content
        i0 = i[2] < 0 ? i[2]+x.parent.rank+1 : i[2]+x.parent.rank
        for j in y.content
            j0 = j[2] < 0 ? j[2]+x.parent.rank+1 : j[2]+x.parent.rank        
            val = x.parent.structure_constants[i0,j0]
            add!( res, [ (i[1]*j[1]*x[1],x[2]) for x in val ]) 
        end
    end
    return lie_algebra_element{Int64}( x.parent, res )
end     

=#

function Base.:*( x::lie_algebra_element, y::lie_algebra_element )::lie_algebra_element

    if x.parent != y.parent 
        throw( "error: not elements of the same algebra" )
    end 

    dom::DataType = domain( x.parent )
    res::Vector{dom} = zeros( dom, dimension( x.parent ))
    xc, yc = x.content, y.content
    sc = x.parent.structure_constants

    dim::Int64 = dimension(x.parent)
    local k::Tuple{dom,Int64}
    local i::Int64, j::Int64
    #local xi::Int64, yi::Int64 

    for i in 1:dim
        xi::Int64, yi::Int64 = xc[i], yc[i]
        for j in i+1:dim
            prod_ij = x.parent.structure_constants[i,j]
            if length( prod_ij ) > 0
                for k in prod_ij
                    res[k[2]] += xi*yc[j]*k[1]
                    res[k[2]] -= yi*xc[j]*k[1]
                end
            end
        end  
    end 

    return lie_algebra_element{domain( x.parent )}( x.parent, res )
end     


function __mult( xc::Vector{Int64}, yc::Vector{Int64}, 
                    sc::Matrix{Vector{Pair{Int64,Int64}}} )::Vector{Int64}

    dim::Int64 = length( xc )                
    res::Vector{Int64} = zeros( Int64, dim )
    local k::Tuple{Int,Int64}
    local i::Int64, j::Int64
    
    for i in 1:dim
        xi::Int64, yi::Int64 = xc[i], yc[i]
        for j in i+1:dim
            prod_ij = sc[i,j]
            if length( prod_ij ) > 0
                for k in prod_ij
                    res[k[2]] += xi*yc[j]*k[1]
                    res[k[2]] -= yi*xc[j]*k[1]
                end
            end
        end  
    end 

    return res
end     

#function  AffineLieAlgebra( type, rank1, rank2, domain )

    
