using Test, ASE, NeighbourLists, JuLIP
using LinearAlgebra: norm
using MolSimPy: neighbourlist

a0 = 2.025

h2("Check neighbourlist without periodicity")
# TODO: implement a test with periodicity as well???
at = bulk("Al", cubic=true) * 3
set_pbc!(at, (false,false,false))
rcut = 1.7 * a0
nlist = neighbourlist(at, rcut)

@info "   ... assemble neighbour list ..."
# create a neighbourlist via a naive double-loop
simple = zeros(length(at), length(at))
X = positions(at)
for n = 2:length(at), m = 1:n-1
   if norm(X[m]-X[n]) <= rcut
      simple[n,m] = simple[m,n] = 1
   end
end

@info "   ... check the bond-iterator ... "
pass_bonds_test = true
for (i,j,r,R) in pairs(nlist)
   if !( (simple[i,j] == 1) && (abs(norm(X[i]-X[j]) - r) < 1e-12) &&
         (norm(X[j]-X[i] - R) < 1e-12)  )
      global pass_bonds_test = false
      break
   end
   # switch the flag
   simple[i,j] = -1
end

h3("check that no false bonds have been found")
println(@test pass_bonds_test)
h3("check that all bonds have been found")
println(@test maximum(simple) == 0)

# revert to original
simple *= -1
h3("   ... check the site iterator ... ")
pass_site_test = true
for (i,j,r,R) in sites(nlist)
   for n = 1:length(j)
      if simple[i,j[n]] != 1
         global pass_site_test = false
         break
      end
      simple[i,j[n]] = -1
   end
   if maximum(simple[i,:]) != 0
      global pass_site_test = false
      break
   end
end
println(@test pass_site_test)
println(@test maximum(simple) == 0)
