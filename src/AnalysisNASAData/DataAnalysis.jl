# Test Traj NASA ###############################################################
# ##############################################################################

# Need further analysis

using CSV
using DataFrames

file = "mars_edl_low_ld_v2.txt"
F = CSV.File(file)
df = DataFrame(F)
names(df)

Plots.plot(df.time, (sqrt.(df.xi.^2+df.yi.^2+df.zi.^2)))

Plots.plot(df.time, (sqrt.(df.vxi.^2+df.vyi.^2+df.vzi.^2)))

Plots.plot(df.time, df.vxi)
Plots.plot(df.time, df.vyi)
Plots.plot(df.time, df.vzi)
Plots.plot(df.time, df.mass) #constant
Plots.plot(df.time, df.piti)
Plots.plot(df.time, df.yawi)
Plots.plot(df.time, df.roli)

Plots.plot(df.time, df.sref) #constant

Plots.plot(df.time, df.pitr)
Plots.plot(df.time, df.rolr)
Plots.plot(df.time, df.yawr)

Plots.plot(df.time, df.alpha) #almost always constant here

Plots.plot(df.time, df.jdate) #constant makes sense
