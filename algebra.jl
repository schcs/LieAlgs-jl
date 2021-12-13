using Base.Threads 
# the type for structure constants tables
sc_type = Matrix{Vector{Tuple{Int64,Int64}}}


# the type for holding information about Lie algebras
struct lie_algebra{T}                     # T is the coefficient domain (Int64, etc)
    type::String                          # type is "A", "B", etc
    rank::Int64                           # Lie rank
    cartan_matrix::Matrix{Int64}          # Cartan matrix
    simple_roots::Vector{Vector{Int64}}   # list of simple roots
    root_system::Vector{Vector{Int64}}    # list of roots
    structure_constants::Vector{Vector{Tuple{Int64,Int64,T}}}        
                                          # the structure constants table 
end

# the type for a Lie algebra element 
mutable struct lie_algebra_element{T}
    parent::lie_algebra                   # the parent Lie algebra
    content::Vector{T}                    # the element
end

# the function for constructing Lie algebras associated with root systems
# and for calculating with their elements
function lie_algebra( type::String, rank::Int64; dom = Int64 )::lie_algebra
    sr::Vector{Vector{Int64}} = simple_root_system( type, rank )
    cm::Matrix{Int64} = cartan_matrix( sr )
    root_system::Vector{Vector{Int64}} = sort( complete_root_system( sr ))
    
    dim::Int64 = rank+length(root_system)
    table::Vector{Vector{Tuple{Int64,Int64,dom}}} = fill( [], dim )
    for i in 1:dim, j in i+1:dim
        ml::Vector{Tuple{dom,Int64}} = mult( sr, root_system,  cm, i, j )
        if length( ml ) > 0
            for m in ml
                push!( table[m[2]], (i,j,m[1]))
            end 
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

function domain( l )
    
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
 
function __mult( xc::Vector{Int64}, yc::Vector{Int64}, 
                    sc::Matrix{Vector{Tuple{Int64,Int64}}} )::Vector{Int64}

    dim::Int64 = length( xc )                
    res::Vector{Int64} = zeros( Int64, dim )
    local k::Tuple{Int,Int64}
    local i::Int64, j::Int64
    xc, yc, sc = x.content, y.content, x.parent.structure_constants
    res = fill( 0*xc[1], length( xc ))
    for (i,j,scij) in sc
          for k in scij
            res[k[2]] += (xc[i]*yc[j]-yc[i]*xc[j])*k[1]
          end  
    end 
    return lie_algebra_element{domain(x.parent)}( x.parent, res )
end  

function __mult__( xc, yc, sc )

    res = fill( 0*xc[1], length( xc ))
    @threads for k2 in 1:length(xc)
        for (i,j,k1) in sc[k2]
            res[k2] += (xc[i]*yc[j]-yc[i]*xc[j])*k1
        end
    end 
    return res
end  


function Base.:*( x::lie_algebra_element, y::lie_algebra_element; 
                sc = x.parent.structure_constants, T = domain(x.parent) )::lie_algebra_element
    
    #= 
    #T::DataType = @domain( x.parent )
    #T = Int64
    xc::Vector{T}, yc::Vector{T} = x.content, y.content
    res::Vector{T} = fill( 0*xc[1], length( xc ))

    for (i,j,scij) in sc
          for k in scij
            res[k[2]] += (xc[i]*yc[j]-yc[i]*xc[j])*k[1]
          end  
    end
    =# 
    res = __mult__( x.content, y.content, sc )
    return lie_algebra_element{T}( x.parent, res )
end  

#=
function mkmult(simpsc)
    function mult(xc, yc)
      dim=length(xc)                
      res=fill(0, dim)
      for (i,j,scij) in simpsc
        for k in scij
          res[k[2]] += (xc[i]*yc[j]-yc[i]*xc[j])*k[1]
        end  
      end 
      res
    end     
  
    mult
  end
  
  using JLD2
  function mkmult(fname::String)
    sc=load(fname)["sc"]
    dim,_=size(sc)
    simpsc=typeof((1,1,sc[1,1]))[]
    for i in 1:dim, j in i+1:dim
      if length(sc[i,j])>0
        push!(simpsc,(i,j,sc[i,j]))
      end
    end
    mkmult(simpsc),dim
  end
  
  mult_,dim=mkmult("sc_table.jld2")
  function comp(K)
    for _ in 1:K
      x,y,z=rand(-10:10,dim),rand(-10:10,dim),rand(-10:10,dim)
      res=mult(mult(x, y),z) + mult(mult(y,z),x) + mult(mult(z, x), y)
      @assert sum(res)==0
    end
  end
  =#