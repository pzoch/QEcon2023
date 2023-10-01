# this is a sample .jl script

using Plots

x = 1:10; y = rand(10); # These are the plotting data
plot(x,y, label="my label")
# save the plot as a .png file
savefig("myplot.png")


