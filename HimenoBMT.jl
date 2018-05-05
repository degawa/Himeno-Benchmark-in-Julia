#This program is Himeno benchmark problem written in Julia 0.6.2.
#This program was simply rewritten into Julia from a program which rewritten in Modern Fortran style.
#https://github.com/degawa/Himeno-Benchmark-in-Modern-Fortran
#
#For the original version of the Himeno benchmark, 
#please refer to the URLs below:
#http://accc.riken.jp/supercom/himenobmt/
#http://accc.riken.jp/en/supercom/himenobmt/
#
#This program is free open-source software distributed under LGPL version 2
#or any later version, inheriting the license of the original version of the Himeno benchmark.

function jacobi!(p, ε, a, b, c, bnd, src, wrk, numItr)
#function jacobi!(p::Array{Float32,3}, ε::Float32, a::Array{Float32,4}, b::Array{Float32,4}, c::Array{Float32,4}, bnd::Array{Float32,3}, src::Array{Float32,3}, wrk::Array{Float32,3}, numItr::Int64)
mimax, mjmax, mkmax = size(p)
    const imax = mimax-1
    const jmax = mjmax-1
    const kmax = mkmax-1

    const ω=0.8f0 #relaxation parameter
    pⁿ⁺¹ = 0.0f0
    Δp   = 0.0f0
    for loop in 1:numItr
        ε = 0.0f0
        for k in 2:kmax-1
        for j in 2:jmax-1
        for i in 2:imax-1
            pⁿ⁺¹ =  a[i,j,k,1 ]*  p[i+1,j  ,k  ]
                   +a[i,j,k,2 ]*  p[i  ,j+1,k  ]
                   +a[i,j,k,3 ]*  p[i  ,j  ,k+1]
                   +b[i,j,k,1 ]*( p[i+1,j+1,k  ]-p[i+1,j-1,k  ] 
                                 -p[i-1,j+1,k  ]+p[i-1,j-1,k  ])
                   +b[i,j,k,2 ]*( p[i  ,j+1,k+1]-p[i  ,j-1,k+1] 
                                 -p[i  ,j+1,k-1]+p[i  ,j-1,k-1])
                   +b[i,j,k,3 ]*( p[i+1,j  ,k+1]-p[i-1,j  ,k+1] 
                                 -p[i+1,j  ,k-1]+p[i-1,j  ,k-1])
                   +c[i,j,k,1 ]*  p[i-1,j  ,k  ]
                   +c[i,j,k,2 ]*  p[i  ,j-1,k  ]
                   +c[i,j,k,3 ]*  p[i  ,j  ,k-1]
                   +src[i,j,k]

            Δp = (pⁿ⁺¹*a[i,j,k,4] - p[i,j,k])*bnd[i,j,k]
            ε += Δp^2
            wrk[i,j,k] = p[i,j,k] + ω*Δp
        end
        end
        end
        p[2:imax-1,2:jmax-1,2:kmax-1] = wrk[2:imax-1,2:jmax-1,2:kmax-1]
    end
end

#Problem size
const XS = (  65,  33,  33)
const S  = ( 129,  65,  65)
const M  = ( 257, 129, 129)
const L  = ( 513, 257, 257)
const XL = (1025, 513, 513)

#Set Grid Parameters
println("Select Grid-size:")
println("           XS (64x32x32)")
println("           S  (128x64x64)")
println("           M  (256x128x128)")
println("           L  (512x256x256)")
println("           XL (1024x512x512)")
print("Grid-size= ")
GridSize = readline(STDIN)          #keyboard input
GridSize = uppercase(GridSize)      #convert to uppercase

if GridSize ≠ "XS" && GridSize ≠ "S" && GridSize ≠ "M" && GridSize ≠ "L" && GridSize ≠ "XL"
    error("Unexpected GridSize $GridSize")
end

#println(eval(Symbol(GridSize)))
GridSize = eval(Symbol(GridSize))   #convert String to Tuple
#println(GridSize)

const mimax = GridSize[1]
const mjmax = GridSize[2]
const mkmax = GridSize[3]

const imax = mimax-1
const jmax = mjmax-1
const kmax = mkmax-1

println("mimax=$mimax mjmax=$mjmax mkmax=$mkmax")
println(" imax=$imax  jmax=$jmax  kmax=$kmax")

#Parameters related to performance measurments
const NanosecToSec = 1e-9
const FlopToMFlop  = 1e-6
const numFlopPerPoint = 34.0
const MFlopsPenIII600 = 82.84

const numPoints = (kmax-2)*(jmax-2)*(imax-2)
const flop = numPoints*numFlopPerPoint

#Declaring and Initializing matrixes
p   = Array{Float32}(mimax,mjmax,mkmax)
a   = Array{Float32}(mimax,mjmax,mkmax,4)
b   = Array{Float32}(mimax,mjmax,mkmax,3)
c   = similar(b)
src = similar(p)
bnd = similar(p)
wrk = similar(p)

a[1:imax,1:jmax,1:kmax,1:3] = 1.0f0
a[1:imax,1:jmax,1:kmax,4  ] = Float32(1.0/6.0)
b[1:imax,1:jmax,1:kmax,:] = 0.0f0
c[1:imax,1:jmax,1:kmax,:] = 1.0f0
bnd[1:imax,1:jmax,1:kmax] = 1.0f0
src[1:imax,1:jmax,1:kmax] = 0.0f0
[ p[:,:,k]=Float32((k-1)^2/(kmax-1)^2) for k = 1:kmax ]

ε = 0.0f0 #error
numItr = 1 #number of iterations of Jacomi method

#warming-up.
#Warm-up is necessary in Julia to corrctly estimate elapsed time excluding JIT compile time
jacobi!(p, ε, a, b, c, bnd, src, wrk, numItr)

#Rehearsal measurment to estimate the number of iterations
numItr = 3
println("Measure the performance in $numItr times.")
println("Start rehearsal measurement process.")

ε = 0.0f0
time_begin_ns = time_ns()
jacobi!(p, ε, a, b, c, bnd, src, wrk, numItr)
time_end_ns   = time_ns()
time_elapsed_s = (time_end_ns-time_begin_ns)*NanosecToSec

mflops = flop*FlopToMFlop / (time_elapsed_s/numItr)
println("  MFLOPS: $mflops  time(s):$time_elapsed_s $ε")
#end Rehearsal measurment

#Acatual measurment
const ExecutionTime_s = 60.0 #[s]
numItr = Int64( div(ExecutionTime_s, time_elapsed_s/numItr) )
println("Now, start the actual measurement process.")
println("The loop will be excuted in $numItr times.")
println("This will take about one minute.")
println("Wait for a while.")

time_begin_ns = time_ns()
jacobi!(p, ε, a, b, c, bnd, src, wrk, numItr)
time_end_ns   = time_ns()
time_elapsed_s = (time_end_ns-time_begin_ns)*NanosecToSec

mflops = flop*FlopToMFlop / (time_elapsed_s/numItr)
score = mflops/MFlopsPenIII600

println(" Loop executed for $numItr times")
println(" Error : $ε")
println(" MFLOPS: $mflops  time(s): $time_elapsed_s")
println(" Score based on Pentium III 600MHz : $score")
#end Acatual measurment
