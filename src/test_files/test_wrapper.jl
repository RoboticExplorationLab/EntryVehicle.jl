#path is the location of the executable
#cd sets the location for working with files
#Note: by setting the cd. The process takes input file in it and write the
#      output file in it as well.
path = "C:\\Users\\33645\\Documents\\Stanford\\AA290ZAC - RA Fall 2019\\MarsGRAM2010-master\\Executables\\marsgram_M10.exe"
cd("C:\\Users\\33645\\Documents\\Stanford\\AA290ZAC - RA Fall 2019\\entry_vehicle\\EntryVehicle.jl\\src\\test_files")

gstdin = Pipe()
gstdout = Pipe()
gstderr = Pipe()
gproc = run(pipeline(`$path`,
                     stdin = gstdin, stdout = gstdout, stderr = gstderr),
            wait = false)
process_running(gproc) || error("There was a problem starting up marsgram.")
w = write(gstdin, string("inputstd6.txt", "\n")) #look for the file in the cd, not in dir of the .exe file
if !(w > 0)
    println("Something went wrong writing to marsgram STDIN.")
    return
end
flush(gstdin)
close(gstdout.in)
close(gstderr.in)
close(gstdin.out)
close(gproc)

#= RESEARCH AND COMMANDS FOR SHELL ACCESS
p = run(`$path`)
run(`inputstd0.txt`)

p = open(pipeline(`$path`), "r+")
write(p, "inputstd6.txt")
close(p)

p = open(path, "r+");
write(p, `inputstd6.txt`); #modifies the .exe or the rights on it for some reasons
close(p)

run(`cmd /c start $path`)

run(pipeline(`$path`))

run(`$path inputstd6.txt`)
=#
