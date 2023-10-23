# this scipt shows you how to work with modules
using Revise # without Revise it is not possible to update the module without restarting the kernel
include("mymodule.jl")




x = 5
y = myfunction(x)

println(y)

