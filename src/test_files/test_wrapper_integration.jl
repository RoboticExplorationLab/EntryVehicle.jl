#Integration of the different steps for MGram wrapper
using DataFrames
using CSV
using ApproxFun

#= PROCESS DECOMPOSITION
write_namelist_file("inputstd8")
execute_marsgram("inputstd8")
data_table = output_data_extraction("OUTPUT2")
=#

function write_namelist_file(name)
    p = open(string(name, ".txt"), "w")
    S = read("namelist_template.txt", String)
    write(p, S)
    close(p)
end

function execute_marsgram(input_file_name)
    path = "C:\\Users\\33645\\Documents\\Stanford\\AA290ZAC - RA Fall 2019\\MarsGRAM2010-master\\Executables\\marsgram_M10.exe"
    cd("C:\\Users\\33645\\Documents\\Stanford\\AA290ZAC - RA Fall 2019\\entry_vehicle\\EntryVehicle.jl\\src\\test_files")

    gstdin = Pipe()
    gstdout = Pipe()
    gstderr = Pipe()
    gproc = run(pipeline(`$path`,
                         stdin = gstdin, stdout = gstdout, stderr = gstderr),
                wait = false)
    process_running(gproc) || error("There was a problem starting up marsgram.")
    w = write(gstdin, string(input_file_name, ".txt\n")) #look for the file in the cd, not in dir of the .exe file
    if !(w > 0)
        println("Something went wrong writing to marsgram STDIN.")
        return
    end
    flush(gstdin)
    close(gstdout.in)
    close(gstderr.in)
    close(gstdin.out)
    close(gproc)
end

function output_data_extraction(output_file_name) #with txt
    data = CSV.File(output_file_name; delim=' ', ignorerepeated=true)
    data_table = DataFrame(data)
    return data_table
end

function write_file(name; LSTFL="INPUT.txt", OUTFL="OUTPUT.txt", PROFILE="null", WAVEFILE="null",
        IERT=1,IUTC=1, MONTH=7, MDAY=20, MYEAR=20, NPOS=41, IHR=12, IMIN=30, SEC=0.0, LONEW=0,
        DUSTTAU=0, DUSTMIN=0.3, DUSTMAX=1.0, DUSTNU=0.003, DUSTDIAM=5.0, DUSTDENS=3000., ALS0=0.0, ALSDUR=48.,
        INTENS=0.0, RADMAX=0.0, DUSTLAT=0.0, DUSTLON=0.0, F107=68.0, NR1=1234, NVARX=1, NVARY=0,
        LOGSCALE=0, FLAT=22.48, FLON=47.97, FHGT=-5., MOLAHGTS=1, HGTASFCM=0.0, ZOFFSET=3.25, IBOUGHER=1,
        DELHGT=0.1, DELLAT=0.5, DELLON=0.5, DELTIME=500.0, DELTATEX=0.0, PROFNEAR=0.0, PROFFAR=0.0, RPSCALE=1.0,
        RWSCALE=1.0, WLSCALE=1.0, WMSCALE=1.0, BLWINFAC=1.0, NMONTE=1, IUP=13, WAVEA0=1.0, WAVEDATE=0.0,
        WAVEA1=0.0, WAVEPHI1=0.0, PHI1DOT=0.0, WAVEA2=0.0, WAVEPHI2=0.0, PHI2DOT=0.0, WAVEA3=0.0,
        WAVEPHI3=0.0, PHI3DOT=0.0, IUWAVE=0, WSCALE=20., CORLMIN=0.0, IPCLAT=1, REQUA=3396.19, RPOLE=3376.20,
        MAPYEAR=0, IDAYDATA=1)

    S =""" \$INPUT_M10\r\n  LSTFL    = '$LSTFL'\r\n  OUTFL    = '$OUTFL'\r\n  profile  = '$PROFILE'\r\n  WaveFile = '$WAVEFILE'\r\n  DATADIR  = 'C:\\Users\\33645\\Documents\\Stanford\\AA290ZAC - RA Fall 2019\\MarsGRAM2010-master\\binFiles\\'\r\n  GCMDIR   = 'C:\\Users\\33645\\Documents\\Stanford\\AA290ZAC - RA Fall 2019\\MarsGRAM2010-master\\binFiles\\'\r\n  IERT     = $IERT\r\n  IUTC     = $IUTC\r\n  MONTH    = $MONTH\r\n  MDAY     = $MDAY\r\n  MYEAR    = $MYEAR\r\n  NPOS     = $NPOS\r\n  IHR      = $IHR\r\n  IMIN     = $IMIN\r\n  SEC      = $SEC\r\n  LonEW    = $LONEW\r\n  Dusttau  = $DUSTTAU\r\n  Dustmin  = $DUSTMIN\r\n  Dustmax  = $DUSTMAX\r\n  Dustnu   = $DUSTNU\r\n  Dustdiam = $DUSTDIAM\r\n  Dustdens = $DUSTDENS\r\n  ALS0     = $ALS0\r\n  ALSDUR   = $ALSDUR\r\n  INTENS   = $INTENS\r\n  RADMAX   = $RADMAX\r\n  DUSTLAT  = $DUSTLAT\r\n  DUSTLON  = $DUSTLON\r\n  MapYear  = $MAPYEAR\r\n  F107     = $F107\r\n  NR1      = $NR1\r\n  NVARX    = $NVARX\r\n  NVARY    = $NVARY\r\n  LOGSCALE = $LOGSCALE\r\n  FLAT     = $FLAT\r\n  FLON     = $FLON\r\n  FHGT     = $FHGT\r\n  MOLAhgts = $MOLAHGTS\r\n  hgtasfcm = $HGTASFCM\r\n  zoffset  = $ZOFFSET\r\n  ibougher = $IBOUGHER\r\n  DELHGT   = $DELHGT\r\n  DELLAT   = $DELLAT\r\n  DELLON   = $DELLON\r\n  DELTIME  = $DELTIME\r\n  deltaTEX = $DELTATEX\r\n  profnear = $PROFNEAR\r\n  proffar  = $PROFFAR\r\n  rpscale  = $RPSCALE\r\n  rwscale  = $RWSCALE\r\n  wlscale  = $WLSCALE\r\n  wmscale  = $WMSCALE\r\n  blwinfac = $BLWINFAC\r\n  NMONTE   = $NMONTE\r\n  iup      = $IUP\r\n  WaveA0   = $WAVEA0\r\n  WaveDate = $WAVEDATE\r\n  WaveA1   = $WAVEA1\r\n  Wavephi1 = $WAVEPHI1\r\n  phi1dot  = $PHI1DOT\r\n  WaveA2   = $WAVEA2\r\n  Wavephi2 = $WAVEPHI2\r\n  phi2dot  = $PHI2DOT\r\n  WaveA3   = $WAVEA3\r\n  Wavephi3 = $WAVEPHI3\r\n  phi3dot  = $PHI3DOT\r\n  iuwave   = $IUWAVE\r\n  Wscale   = $WSCALE\r\n  corlmin  = $CORLMIN\r\n  ipclat   = $IPCLAT\r\n  requa    = $REQUA\r\n  rpole    = $RPOLE\r\n  idaydata = $IDAYDATA\r\n \$END\r\n\r\n """
    p = open(string(name, ".txt"), "w")
    write(p, S)
    close(p)
end

function marsgram_wrapper(name; kwargs...)
    write_file(name; kwargs...)
    execute_marsgram(name)
    sleep(1.0) #seems to be necessary otherwise, no time to create the file. This is a problem here
    data_table = output_data_extraction(kwargs[:OUTFL]) #okay that works because kwargs is a dictionnary actually. Awesome
    return data_table
end

X = marsgram_wrapper("test9"; OUTFL = "OUT21.txt", NPOS=100)
