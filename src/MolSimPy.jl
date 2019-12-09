module MolSimPy

export MatSciPy

include("MatSciPy.jl")

using .MatSciPy
using ASE: ASEAtoms

function neighbourlist(at::ASEAtoms, cutoff::Float64)::MatSciPy.NeighbourList
   # TODO: also recompute if rcut is different !!!!!
   # if no previous neighbourlist is available, compute a new one
   # if !has_transient(at, (:nlist, cutoff)) || recompute
   #    # this nlist will be destroyed as soon as positions change
   #    set_transient!(at, (:nlist, cutoff), MatSciPy.NeighbourList(at, cutoff))
   # end
   # return get_transient(at, (:nlist, cutoff))
   return MatSciPy.NeighbourList(at, cutoff)
end

end # module
