using DataFrames
using CSV
using Plots
pyplot()

X = CSV.read(joinpath(pwd(), "OUTPUT.txt"))
show(X)
describe(X)
dropmissing(X)

XX = CSV.File(joinpath(pwd(), "OUTPUT.txt"))
XXX = CSV.File(joinpath(pwd(),"OUTPUT.txt"); delim=' ', ignorerepeated=true)

show(XXX)
XXX = DataFrame(XXX)
XXX.Time
XXX.Denkgm3

#Density with respect to Altitude for the specified conditions
plot(XXX.Denkgm3, XXX.HgtMOLA)
