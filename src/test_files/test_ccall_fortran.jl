#test ccall Julia. C/Fortran interfacing
using Libdl

ccall((:MarsGRAM_M10, "marsgram"), Int32, ())
cd("C:\\Users\\33645\\Documents\\Stanford\\AA290ZAC - RA Fall 2019\\MarsGRAM2010-master\\Code")
push!(Libdl.DL_LOAD_PATH,"C:\\Users\\33645\\Documents\\Stanford\\AA290ZAC - RA Fall 2019\\MarsGRAM2010-master\\Code\\marsgram.so")
path = pwd()
file = joinpath("$path", "inputstd6.txt")

typeof(file)

#test
#library should be libc but actually here for some reasons on windows msvcrt works better
t = ccall((:clock, "msvcrt"), Int32, ()) #need the libc.so file and need to work in the same directory.

ccall((:GetDoubleClickTime, "User32"), stdcall, Int32, (), )

pwd()
ccall((:atmos2_m10_, "marsgram_test"), Int32, ())


path = "/home/rexlab/remy/MarsGRAM2010-master/Code"
cd("$path")
Libdl.dlopen("/home/rexlab/remy/MarsGRAM2010-master/Code/marsgram_test.so")
ccall((:atmos2_m10_, "/home/rexlab/remy/MarsGRAM2010-master/Code/marsgram_test.so"), Int32, ())

pwd()


A = 1.0
B = 1.0
C = 1.0
D = 2.0
E = 3.0
F = 1.0
G = 1.0
H = 1.0
I = 0.5
J = 5.0
Ref{A}
Ptr{A}
Ref(Float64)

path = "/home/rexlab/remy/MarsGRAM2010-master/Code"
cd("$path")
marsgram = Libdl.dlopen("/home/rexlab/remy/MarsGRAM2010-master/Code/marsgram_test.so")
ccall(Libdl.dlsym(marsgram,:bltp_m10_),
    Nothing, (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
    Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}), A, B, C, D, E, F, G, H, I, J)
println(J)

#testing with my own fortran function here

ftest = Libdl.dlopen("/home/rexlab/remy/MarsGRAM2010-master/Code/test.so") #recognize my library
t = ones(Int32,2)
ccall((:foo1_, "/home/rexlab/remy/MarsGRAM2010-master/Code/test.so"), Nothing, (Ref{Int32}, Ref{Int32}, Ptr{Int32}), Int32(2), Int32.([3, -1]), t)
println(t) #Nice it works so I can call fortran code now, that's neat, even with arrays and stuff. Now what about the code I have, MARSGRAM.

t = ones(Int64,2)
ccall(Libdl.dlsym(ftest,:foo2_), Nothing, (Ref{Int64}, Ref{Int64}, Ptr{Int64}), 2, [3, -1], t) #okay nice, that's another way to call my thing here and it works pretty good (in the library), nice
ccall((:foo2_, "/home/rexlab/remy/MarsGRAM2010-master/Code/test.so"), Nothing, (Ref{Int64}, Ref{Int64}, Ptr{Int64}), 2, [3, -1], t)
println(t)

#okay functions should be before subroutines in the f90 code actually.
ftest = Libdl.dlopen("/home/rexlab/remy/MarsGRAM2010-master/Code/test.so")
x1 = 7
a1 = [0]
b1 = [0]
r1 = ccall(Libdl.dlsym(ftest, :foo_), Int64,
            (Ref{Int64},), x1)
println(r1)

x2 = 7.0
a2 = Cdouble[1.0]
b2 = Cdouble[1.0]
ccall(Libdl.dlsym(ftest, :keg_), Nothing,
      (Ref{Float64}, Ref{Float64}, Ref{Float64}), x2, a2, b2)

println(a2[1])
println(b2[1])

#nice I understand, for subroutine you basically return Nothing all the time,
#as you only modify stuff in inputs etc...



####################################
######### test MARSGRAM 1 ##########
######### checking calls ###########
######### using .txt input##########
####################################


#test with MARSGRAM_M10 functions/subroutines:

#SUBROUTINE SETUP_M10(CHGT, CLAT, CLON, CSEC, DAY0, RHOD, RHOU, RHOV, RHOW
#   , DHGT, DLAT, DLON, DTIME, MAXNUM, NRN1, NMCR1, DSUNLS, DRADAU, DOWLT
#   , LNEW, INPUTFL, IUSTDOUT, IULIST, HGTASFC, INERT, INUTC, STEPMIN,
#   PROFNR, PROFFR, NPROF)



#Input data to MarsGRAM
LSTFL="INPUT.txt"; OUTFL="OUTPUT.txt"; PROFILE="null"; WAVEFILE="null";
IERT=1;IUTC=1; MONTH=7; MDAY=20; MYEAR=20; NPOS=41; IHR=12; IMIN=30; SEC=0.0; LONEW=0;
DUSTTAU=0; DUSTMIN=[0.3]; DUSTMAX=[1.0]; DUSTNU=[0.003]; DUSTDIAM=[5.0]; DUSTDENS=3000.; ALS0=[0.0]; ALSDUR=[48.];
INTENS=[0.0]; RADMAX=[0.0]; DUSTLAT=[0.0]; DUSTLON=[0.0]; F107=[68.0]; NR1=1234; NVARX=1; NVARY=0;
LOGSCALE=0; FLAT=[22.48]; FLON=[47.97]; FHGT=[-5.]; MOLAHGTS=1; HGTASFCM=[0.0]; ZOFFSET=[3.25]; IBOUGHER=1;
DELHGT=[0.1]; DELLAT=[0.5]; DELLON=[0.5]; DELTIME=[500.0]; DELTATEX=0.0; PROFNEAR=0.0; PROFFAR=0.0; RPSCALE=1.0;
RWSCALE=[1.0]; WLSCALE=[1.0]; WMSCALE=[1.0]; BLWINFAC=[1.0]; NMONTE=1; IUP=13; WAVEA0=[1.0]; WAVEDATE=[0.0];
WAVEA1=[0.0]; WAVEPHI1=[0.0]; PHI1DOT=[0.0]; WAVEA2=[0.0]; WAVEPHI2=[0.0]; PHI2DOT=[0.0]; WAVEA3=[0.0];
WAVEPHI3=[0.0]; PHI3DOT=[0.0]; IUWAVE=0; WSCALE=[20.]; CORLMIN=[0.0]; IPCLAT=1; REQUA=[3396.19]; RPOLE=[3376.20];
MAPYEAR=0 ;IDAYDATA=1


CHGT = [1.0]; CLAT = [1.0]; CLON = [1.0]; CSEC = [1.0]; DAY0 = [1.0]; RHOD=[1.0]; RHOU=[1.0]; RHOV=[1.0]; RHOW=[1.0]
     ; DHGT=[1.0]; DLAT=[1.0]; DLON=[1.0]; DTIME=[1.0]; MAXNUM = [1]; NRN1=[1]; NMCR1=[1]; DSUNLS=[1.0]; DRADAU=[1.0]; DOWLT=[1.0]
     ; LNEW=[1]; INPUTFL = "inputstd6.txt"; IUSTDOUT = 6; IULIST=[1]; HGTASFC=[1.0]; INERT=[1]; INUTC=[1]; STEPMIN=[1.0];
     PROFNR=[1.0]; PROFFR=[1.0]; NPROF=[1]

using Libdl
push!(Libdl.DL_LOAD_PATH,"/home/rexlab/remy/MarsGRAM2010-master/Code/marsgram_test.so")

path = "/home/rexlab/remy/MarsGRAM2010-master/Code"
cd("$path")
marsgram = Libdl.dlopen("/home/rexlab/remy/MarsGRAM2010-master/Code/marsgram_test.so") #open shared library
ccall(Libdl.dlsym(marsgram,:setup_m10_),
    Nothing, (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
    Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
    Ref{Float64}, Ref{Float64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Float64},
    Ref{Float64}, Ref{Int64}, Ptr{UInt8}, Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Int64}, Ref{Int64},
    Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64}), CHGT, CLAT, CLON, CSEC, DAY0, RHOD, RHOU, RHOV, RHOW,
    DHGT, DLAT, DLON, DTIME, MAXNUM, NRN1, NMCR1, DSUNLS, DRADAU, DOWLT, LNEW, INPUTFL, IUSTDOUT, IULIST, HGTASFC,
    INERT, INUTC, STEPMIN, PROFNR, PROFFR, NPROF)

#Save initial position, date, and position displacement values

#initialize
NUMWAVE = 0;
PERTSTEP =[0.0];
IUPDATE = 0;
EOF = 0;

# ALl these things are OUT DOUBLE (so simply give them a value here)
TEMP =[1.0]; PRES =[1.0]; DENSLO=[1.0];DENS=[1.0];DENSHI=[1.0]; DENSP=[1.0]; EWWIND=[1.0]; EWPERT=[1.0]; NSWIND=[1.0]; NSPERT=[1.0];
VWPERT=[1.0]; HRHO=[1.0]; HPRES=[1.0];
CORLIM=[1.0]; DENSTOT=[1.0]; ALS=[1.0]; SZA=[1.0]; OWLT=[1.0]; SUNLAT=[1.0]; SUNLON=[1.0]; MARSAU=[1.0]; TLOCAL=[1.0];

#Initialize density array
dens_results = []

#Modifications expected if MonteCarlo Number is more than 1
for I=0:1:MAXNUM[1]
    ccall(Libdl.dlsym(marsgram,:datastep_m10_), Nothing, (Ref{Int64}, Ref{Float64},
    Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64},
    Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
    Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
    Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
    Ref{Float64}, Ref{Float64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Int64}, Ref{Float64}, Ref{Int64},
    Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
    Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64}), I, CHGT, CLAT, CLON, CSEC, DAY0,
    RHOD, RHOU, RHOV, RHOW, EOF, DELHGT, DELLAT, DELLON, DELTIME, TEMP, PRES,
    DENSLO, DENS, DENSHI, DENSP, EWWIND, EWPERT, NSWIND, NSPERT, VWPERT, HRHO, HPRES,
    0.0, 0.0, 0.0, 0.0, 0.0, LONEW, CORLIM, DENSTOT, NUMWAVE, HGTASFC, IERT, IUTC, PERTSTEP,
    CORLMIN, IUPDATE, ALS, SZA, OWLT, SUNLAT, SUNLON, MARSAU, TLOCAL, PROFNEAR, PROFFAR, NPROF)
    append!(dens_results, DENS)
end



###################
######TEST 2#######
###################
#Removed Namelist and writing files stuff
#removed the INPUTFILE argument in setup_m10_
#ONLY THING REMAINING IS THE DATADIR STUFF THAT IS NOT ACCEPTED BY SETUP (NOt same type of string/character)
#pretty sure it is a problem of type (character or else)

using Libdl
push!(Libdl.DL_LOAD_PATH,"/home/rexlab/remy/MarsGRAM2010-master/Code/marsgram2.so")
path = "/home/rexlab/remy/MarsGRAM2010-master/Code"
cd("$path")
marsgram2 = Libdl.dlopen("/home/rexlab/remy/MarsGRAM2010-master/Code/marsgram2.so") #open shared library

#Former Namelist file arguments
LSTFL = "INPUT.txt"; OUTFL = "OUTPUT.txt"; PROFILE="null"; WAVEFILE="null"; DATADIR="/home/rexlab/remy/MarsGRAM2010-master/binFiles/";
GCMDIR ="/home/rexlab/remy/MarsGRAM2010-master/binFiles/";
IERT=1;IUTC=1; MONTH=7; MDAY=20; MYEAR=20; NPOS=41; IHR=12; IMIN=30; SEC=0.0; LONEW=0;
DUSTTAU=0; DUSTMIN=[0.3]; DUSTMAX=[1.0]; DUSTNU=[0.003]; DUSTDIAM=[5.0]; DUSTDENS=3000.; ALS0=[0.0]; ALSDUR=[48.];
INTENS=[0.0]; RADMAX=[0.0]; DUSTLAT=[0.0]; DUSTLON=[0.0]; F107=[68.0]; NR1=1234; NVARX=1; NVARY=0;
LOGSCALE=0; FLAT=[22.48]; FLON=[47.97]; FHGT=[-5.]; MOLAHGTS=1; HGTASFCM=[0.0]; ZOFFSET=[3.25]; IBOUGHER=1;
DELHGT=[0.1]; DELLAT=[0.5]; DELLON=[0.5]; DELTIME=[500.0]; DELTATEX=0.0; PROFNEAR=0.0; PROFFAR=0.0; RPSCALE=1.0;
RWSCALE=[1.0]; WLSCALE=[1.0]; WMSCALE=[1.0]; BLWINFAC=[1.0]; NMONTE=1; IUP=13; WAVEA0=[1.0]; WAVEDATE=[0.0];
WAVEA1=[0.0]; WAVEPHI1=[0.0]; PHI1DOT=[0.0]; WAVEA2=[0.0]; WAVEPHI2=[0.0]; PHI2DOT=[0.0]; WAVEA3=[0.0];
WAVEPHI3=[0.0]; PHI3DOT=[0.0]; IUWAVE=0; WSCALE=[20.]; CORLMIN=[0.0]; IPCLAT=1; REQUA=[3396.19]; RPOLE=[3376.20];
MAPYEAR=0 ;IDAYDATA=1

CHGT = [1.0]; CLAT = [1.0]; CLON = [1.0]; CSEC = [1.0]; DAY0 = [1.0]; RHOD=[1.0]; RHOU=[1.0]; RHOV=[1.0]; RHOW=[1.0]
     ; DHGT=[1.0]; DLAT=[1.0]; DLON=[1.0]; DTIME=[1.0]; MAXNUM = [1]; NRN1=[1]; NMCR1=[1]; DSUNLS=[1.0]; DRADAU=[1.0]; DOWLT=[1.0]
     ; LNEW=[1]; IUSTDOUT = 6; IULIST=[1]; HGTASFC=[1.0]; INERT=[1]; INUTC=[1]; STEPMIN=[1.0];
     PROFNR=[1.0]; PROFFR=[1.0]; NPROF=[1]

ccall(Libdl.dlsym(marsgram2,:setup_m10_),
    Nothing, (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
    Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
    Ref{Float64}, Ref{Float64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Float64},
    Ref{Float64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Int64}, Ref{Int64},
    Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64}, #= =#Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8},
    Ref{Int64}, Ref{Int64},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Float64},Ref{Int64},Ref{Int64}, Ref{Float64},Ref{Float64},
    Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},
    Ref{Float64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64},
    Ref{Float64},Ref{Float64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},Ref{Float64}, Ref{Float64},
    Ref{Float64}, Ref{Float64},Ref{Float64}, Ref{Float64},Ref{Float64}, Ref{Float64},
    Ref{Int64}, Ref{Int64}, Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},
    Ref{Float64},Ref{Float64},Ref{Float64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Int64}, Ref{Float64}, Ref{Float64},
    Ref{Int64}, Ref{Int64}), CHGT, CLAT, CLON, CSEC, DAY0, RHOD, RHOU, RHOV, RHOW,
    DHGT, DLAT, DLON, DTIME, MAXNUM, NRN1, NMCR1, DSUNLS, DRADAU, DOWLT, LNEW, IUSTDOUT, IULIST, HGTASFC,
    INERT, INUTC, STEPMIN, PROFNR, PROFFR, NPROF, LSTFL, OUTFL, PROFILE, WAVEFILE, DATADIR, GCMDIR, IERT, IUTC, MONTH, MDAY, MYEAR, NPOS
    ,IHR, IMIN, SEC, LONEW, DUSTTAU, DUSTMIN, DUSTMAX, DUSTNU, DUSTDIAM, DUSTDENS
    ,ALS0, ALSDUR, INTENS, RADMAX, DUSTLAT, DUSTLON, F107, NR1, NVARX, NVARY, LOGSCALE
    ,FLAT, FLON, FHGT, MOLAHGTS, HGTASFCM, ZOFFSET, IBOUGHER, DELHGT, DELLAT, DELLON, DELTIME
    , DELTATEX, PROFNEAR, PROFFAR, RPSCALE, RWSCALE, WLSCALE, WMSCALE, BLWINFAC, NMONTE, IUP
    , WAVEA0, WAVEDATE, WAVEA1, WAVEPHI1, PHI1DOT, WAVEA2, WAVEPHI2, PHI2DOT, WAVEA3, WAVEPHI3
    ,PHI3DOT, IUWAVE, WSCALE, CORLMIN, IPCLAT, REQUA, RPOLE, MAPYEAR, IDAYDATA)

#Save initial position, date, and position displacement values

#initialize
NUMWAVE = 0;
PERTSTEP =[0.0];
IUPDATE = 0;
EOF = 0;

# ALl these things are OUT DOUBLE (so simply give them a value here)
TEMP =[1.0]; PRES =[1.0]; DENSLO=[1.0];DENS=[1.0];DENSHI=[1.0]; DENSP=[1.0]; EWWIND=[1.0]; EWPERT=[1.0]; NSWIND=[1.0]; NSPERT=[1.0];
VWPERT=[1.0]; HRHO=[1.0]; HPRES=[1.0];
CORLIM=[1.0]; DENSTOT=[1.0]; ALS=[1.0]; SZA=[1.0]; OWLT=[1.0]; SUNLAT=[1.0]; SUNLON=[1.0]; MARSAU=[1.0]; TLOCAL=[1.0];

#Initialize density array
dens_results = []

#Modifications expected if MonteCarlo Number is more than 1
for I=0:1:MAXNUM[1]
    ccall(Libdl.dlsym(marsgram2,:datastep_m10_), Nothing, (Ref{Int64}, Ref{Float64},
    Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64},
    Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
    Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
    Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
    Ref{Float64}, Ref{Float64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Int64}, Ref{Float64}, Ref{Int64},
    Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
    Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64}), I, CHGT, CLAT, CLON, CSEC, DAY0,
    RHOD, RHOU, RHOV, RHOW, EOF, DELHGT, DELLAT, DELLON, DELTIME, TEMP, PRES,
    DENSLO, DENS, DENSHI, DENSP, EWWIND, EWPERT, NSWIND, NSPERT, VWPERT, HRHO, HPRES,
    0.0, 0.0, 0.0, 0.0, 0.0, LONEW, CORLIM, DENSTOT, NUMWAVE, HGTASFC, IERT, IUTC, PERTSTEP,
    CORLMIN, IUPDATE, ALS, SZA, OWLT, SUNLAT, SUNLON, MARSAU, TLOCAL, PROFNEAR, PROFFAR, NPROF)
    append!(dens_results, DENS)
end

using Libdl
marsgram2 = Libdl.dlopen("/home/rexlab/remy/MarsGRAM2010-master/Code/marsgram2.so")

function wrapper_density(LSTFL="INPUT.txt", OUTFL="OUTPUT.txt",
        PROFILE="null", WAVEFILE="null", DATADIR ="null", GCMDIR="null",
        IERT=1,IUTC=1, MONTH=7, MDAY=20, MYEAR=20, NPOS=41, IHR=12, IMIN=30,
        SEC=0.0, LONEW=0, DUSTTAU=0, DUSTMIN=0.3, DUSTMAX=1.0, DUSTNU=0.003,
        DUSTDIAM=5.0, DUSTDENS=3000., ALS0=0.0, ALSDUR=48., INTENS=0.0,
        RADMAX=0.0, DUSTLAT=0.0, DUSTLON=0.0, F107=68.0, NR1=1234, NVARX=1, NVARY=0,
        LOGSCALE=0, FLAT=22.48, FLON=47.97, FHGT=-5., MOLAHGTS=1, HGTASFCM=0.0,
        ZOFFSET=3.25, IBOUGHER=1, DELHGT=0.1, DELLAT=0.5, DELLON=0.5, DELTIME=500.0,
        DELTATEX=0.0, PROFNEAR=0.0, PROFFAR=0.0, RPSCALE=1.0, RWSCALE=1.0, WLSCALE=1.0,
        WMSCALE=1.0, BLWINFAC=1.0, NMONTE=1, IUP=13, WAVEA0=1.0, WAVEDATE=0.0, WAVEA1=0.0,
        WAVEPHI1=0.0, PHI1DOT=0.0, WAVEA2=0.0, WAVEPHI2=0.0, PHI2DOT=0.0, WAVEA3=0.0,
        WAVEPHI3=0.0, PHI3DOT=0.0, IUWAVE=0, WSCALE=20., CORLMIN=0.0, IPCLAT=1, REQUA=3396.19,
        RPOLE=3376.20, MAPYEAR=0, IDAYDATA=1)

        #Variable Initialization for setup subroutine
        CHGT = [1.0]; CLAT = [1.0]; CLON = [1.0]; CSEC = [1.0]; DAY0 = [1.0]; RHOD=[1.0];
        RHOU=[1.0]; RHOV=[1.0]; RHOW=[1.0]; DHGT=[1.0]; DLAT=[1.0]; DLON=[1.0]; DTIME=[1.0];
        MAXNUM = [1]; NRN1=[1]; NMCR1=[1]; DSUNLS=[1.0]; DRADAU=[1.0]; DOWLT=[1.0];
        LNEW=[1]; IUSTDOUT = 6; IULIST=[1]; HGTASFC=[1.0]; INERT=[1]; INUTC=[1]; STEPMIN=[1.0];
        PROFNR=[1.0]; PROFFR=[1.0]; NPROF=[1]

        #call setup fotran subroutine, setup variables for the run
        ccall(Libdl.dlsym(marsgram2,:setup_m10_),
            Nothing, (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
            Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
            Ref{Float64}, Ref{Float64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Float64},
            Ref{Float64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Int64}, Ref{Int64},
            Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64}, #= =#Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8},
            Ref{Int64}, Ref{Int64},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Float64},Ref{Int64},Ref{Int64}, Ref{Float64},Ref{Float64},
            Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},
            Ref{Float64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64},
            Ref{Float64},Ref{Float64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},Ref{Float64}, Ref{Float64},
            Ref{Float64}, Ref{Float64},Ref{Float64}, Ref{Float64},Ref{Float64}, Ref{Float64},
            Ref{Int64}, Ref{Int64}, Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},
            Ref{Float64},Ref{Float64},Ref{Float64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Int64}, Ref{Float64}, Ref{Float64},
            Ref{Int64}, Ref{Int64}), CHGT, CLAT, CLON, CSEC, DAY0, RHOD, RHOU, RHOV, RHOW,
            DHGT, DLAT, DLON, DTIME, MAXNUM, NRN1, NMCR1, DSUNLS, DRADAU, DOWLT, LNEW, IUSTDOUT, IULIST, HGTASFC,
            INERT, INUTC, STEPMIN, PROFNR, PROFFR, NPROF, LSTFL, OUTFL, PROFILE, WAVEFILE, DATADIR, GCMDIR, IERT, IUTC, MONTH, MDAY, MYEAR, NPOS
            ,IHR, IMIN, SEC, LONEW, DUSTTAU, DUSTMIN, DUSTMAX, DUSTNU, DUSTDIAM, DUSTDENS
            ,ALS0, ALSDUR, INTENS, RADMAX, DUSTLAT, DUSTLON, F107, NR1, NVARX, NVARY, LOGSCALE
            ,FLAT, FLON, FHGT, MOLAHGTS, HGTASFCM, ZOFFSET, IBOUGHER, DELHGT, DELLAT, DELLON, DELTIME
            , DELTATEX, PROFNEAR, PROFFAR, RPSCALE, RWSCALE, WLSCALE, WMSCALE, BLWINFAC, NMONTE, IUP
            , WAVEA0, WAVEDATE, WAVEA1, WAVEPHI1, PHI1DOT, WAVEA2, WAVEPHI2, PHI2DOT, WAVEA3, WAVEPHI3
            ,PHI3DOT, IUWAVE, WSCALE, CORLMIN, IPCLAT, REQUA, RPOLE, MAPYEAR, IDAYDATA)

        #initialize values for
        NUMWAVE = 0; PERTSTEP =[0.0]; IUPDATE = 0; EOF = 0; TEMP =[1.0]; PRES =[1.0];
        DENSLO=[1.0];DENS=[1.0];DENSHI=[1.0]; DENSP=[1.0]; EWWIND=[1.0]; EWPERT=[1.0];
        NSWIND=[1.0]; NSPERT=[1.0]; VWPERT=[1.0]; HRHO=[1.0]; HPRES=[1.0]; CORLIM=[1.0];
        DENSTOT=[1.0]; ALS=[1.0]; SZA=[1.0]; OWLT=[1.0]; SUNLAT=[1.0]; SUNLON=[1.0];
        MARSAU=[1.0]; TLOCAL=[1.0];

        #Initialize density array
        dens_results = []
        height = FHGT:DELHGT:FHGT+(NPOS-1)*DELHGT

        #Modifications expected if MonteCarlo Number is more than 1
        for I=0:1:MAXNUM[1]
            ccall(Libdl.dlsym(marsgram2,:datastep_m10_), Nothing, (Ref{Int64}, Ref{Float64},
            Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64},
            Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
            Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
            Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
            Ref{Float64}, Ref{Float64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Int64}, Ref{Float64}, Ref{Int64},
            Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
            Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64}), I, CHGT, CLAT, CLON, CSEC, DAY0,
            RHOD, RHOU, RHOV, RHOW, EOF, DELHGT, DELLAT, DELLON, DELTIME, TEMP, PRES,
            DENSLO, DENS, DENSHI, DENSP, EWWIND, EWPERT, NSWIND, NSPERT, VWPERT, HRHO, HPRES,
            0.0, 0.0, 0.0, 0.0, 0.0, LONEW, CORLIM, DENSTOT, NUMWAVE, HGTASFC, IERT, IUTC, PERTSTEP,
            CORLMIN, IUPDATE, ALS, SZA, OWLT, SUNLAT, SUNLON, MARSAU, TLOCAL, PROFNEAR, PROFFAR, NPROF)
            append!(dens_results, DENS)
        end

        return height, dens_results
end

h, d = wrapper_density()
#remove writing when using setup_m10_ subroutine
