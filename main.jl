#Load Packages to Be Used
using Printf
using Plots
using DelimitedFiles
using LinearAlgebra
using Dates
using TimerOutputs
using JSON3
using BenchmarkTools

const to = TimerOutput()
@timeit to "Setup" begin 

    include("indicies_scalar.jl")
    include("indicies_field.jl")
    include("equilibrium_o.jl")
    include("collide_CM_thermal.jl")
    include("stream.jl")
    include("cavityint_d2q9.jl")
    include("macrovar.jl")
    include("vout.jl")
    include("writemacro.jl")

    cd(@__DIR__)
    #D2Q9 lattice setup
    ndir = 9
    w = 4.0/9.0
    ws = 1.0 / 9.0
    wd = 1.0 / 36.0
    wi = [w,ws,ws,ws,ws,wd,wd,wd,wd]
    #Domain and array size setup
    lnx = 256
    lny = 256
    nnx = lnx 
    nny = lny 
    nnode = nnx * nny
    popnode = nnode * ndir
    v0 = 1.0
    alx = lnx -1
    aly = lny -1
    dx = alx / (nnx-1)
    dy = alx / (nny-1)
    dt = dx / v0

    #Lattice Speed of sound definition, and combination optimizations
    cs = sqrt(1.0 / 3.0) * v0
    cs2 = cs * cs
    invcs2 = 1.0 / cs2
    invcs2_2 = 0.5 * invcs2
    invcs4 = invcs2 * invcs2
    invcs4_2 = 0.5 * invcs4

    #Paramater Definitions
    #Mach Number
    Ma = 0.01
    #Prandtl Number
    Pr = 0.71
    #Rayleigh Number 
    Ra = 10^3
    #Initial Density
    rhoo = 1.0 

    #Lattice - Analog Calculation
    LLP = dx/lnx
    Lo = alx
    Uo = cs * Ma
    To = Lo / Uo
    Fo = (Uo * Uo) / Lo
    #Set Wall Temperatures
    thot = 2.0 
    tcold = 1.0 
    tmid = 0.5 * (thot + tcold)
    dT = thot - tcold

    #Coeffecent of Thermal Expansion Calculation
    betah = (Uo * Uo) / (Fo * (thot - tcold) * Lo)

    #Thermal Diffusivity Calculation
    alphah = sqrt((Fo * betah * (thot - tcold) * Lo * Lo * Lo)/(Ra * Pr))
    #Kinematic Viscosity
    nu = Pr * alphah
    #Dynamic Viscosity
    mu = rhoo * nu 
    #Reynolds Number
    Re = rhoo * Uo * Lo / mu

    #Relaxation time calculations
    tauf = invcs2 * nu + 0.5
    taut = invcs2 * alphah + 0.5
    #Optimizations of relaxation time calculations
    omegat = 1.0 / taut
    w3 = 1.0
    w4 = 1.0 / tauf 
    phi = w3/w4
    phip = phi + 1.0
    phim = phi - 1.0
    mphi = 1.0 - phi
    invphip = 1.0 / phip
    xo = 0.5 * (lnx-1)

    #Wall Velocities and non-gravity body forces
    uwxr = 0.0
    uwxl = 0.0
    uwxb = 0.0
    uwxt = 0.0
    uwyr = 0.0
    uwyl = 0.0
    uwyb = 0.0
    uwyt = 0.0
    fbx = 0.0
    fby = 0.0

    #Error checking for loop exit conditions. 
    global errorx = 1.0
    global errory = 1.0
    global errort = 1.0

    global xconv = 1
    global yconv = 1
    global tconv = 1

    #Maximum number of iterations allowed
    global nsteps=3000000

    #Printing Paramaters befor simulations
    println("dx= ",dx)
    println("dt= ",dt)
    println("Re= ",Re)
    println("Ra= ",Ra)
    println("nu= ",nu)
    println("mu= ",mu)
    println("Ma= ",Ma)
    println("beta= ",betah)
    println("alpha= ",alphah)
    println("u0= ",Uo)
    println("tauf= ",tauf)
    println("taut= ",taut)
    println("To= ",To)
    println("Fo= ",Fo)
    println("Lo= ",Lo)

    #initilatation of Storage Vectors
    f1 = Vector{Float64}(undef,popnode)
    t1 = Vector{Float64}(undef,popnode)
    f2 = Vector{Float64}(undef,popnode)
    t2 = Vector{Float64}(undef,popnode)
    fcopy = Vector{Float64}(undef,popnode)
    tcopy = Vector{Float64}(undef,popnode)
    solidnode = Vector{Int64}(undef,nnode)
    solidstr = Vector{Int64}(undef,nnode)
    pres = Vector{Float64}(undef,nnode)
    temp = Vector{Float64}(undef,nnode)
    rho = Vector{Float64}(undef,nnode)
    psi = Vector{Float64}(undef,nnode)
    rhoold = Vector{Float64}(undef,nnode)
    tempold = Vector{Float64}(undef,nnode)
    ui = Vector{Float64}(undef,nnode)
    uj = Vector{Float64}(undef,nnode)
    uiold = Vector{Float64}(undef,nnode)
    ujold = Vector{Float64}(undef,nnode)
    fri = Vector{Float64}(undef,nnode)
    frj = Vector{Float64}(undef,nnode)

    #initilatation of Storage Arrays
    printui = Array{Float64}(undef,nnx,nny)
    printuj = Array{Float64}(undef,nnx,nny)
    printrho = Array{Float64}(undef,nnx,nny)
    printtemp = Array{Float64}(undef,nnx,nny)
    solidnodea = Array{Float64}(undef,nnx,nny)
    vp = Array{Float64}(undef,ndir,2)
    vpi = Array{Int64}(undef,ndir,2)

    #Discrete particle velocity directions in both integer and float forms; 1=x 2=y
    vp[1,1] = 0.0     
    vp[2,1] = 1.0     
    vp[3,1] = 0.0     
    vp[4,1] = -1.0    
    vp[5,1] = 0.0    
    vp[6,1] = 1.0    
    vp[7,1] = -1.0    
    vp[8,1] = -1.0    
    vp[9,1] = 1.0     

    vp[1,2] = 0.0
    vp[2,2] = 0.0
    vp[3,2] = 1.0
    vp[4,2] = 0.0
    vp[5,2] = -1.0
    vp[6,2] = 1.0
    vp[7,2] = 1.0 
    vp[8,2] = -1.0
    vp[9,2] = -1.0

    vpi[1,1] = 0     
    vpi[2,1] = 1     
    vpi[3,1] = 0     
    vpi[4,1] = -1    
    vpi[5,1] = 0    
    vpi[6,1] = 1    
    vpi[7,1] = -1    
    vpi[8,1] = -1    
    vpi[9,1] = 1     

    vpi[1,2] = 0
    vpi[2,2] = 0
    vpi[3,2] = 1
    vpi[4,2] = 0
    vpi[5,2] = -1
    vpi[6,2] = 1
    vpi[7,2] = 1 
    vpi[8,2] = -1
    vpi[9,2] = -1

    #Initilazation of Macrovariable feilds, and establishment of domain boundaries; nodes which are boundaries are denoted
    cavityint_d2q9(nnx,nny,rhoo,rho,rhoold,ui,uiold,uj,ujold,pres,temp,tempold,tmid,fri,frj,solidnode,solidstr,cs2)

    #Population initilatation based on Inital Macrovariable feilds
    equilibrium_o(nnx,nny,nnode,ui,uj,temp,rho,ndir,vp,wi,invcs2,invcs4_2,invcs2_2,f1,t1,f2,t2)

    #Creation of output directory for storage of files. 
    time = Dates.format(now(), "YYYYmmdd-HHMM")
    name=joinpath("data_$time RA_$Ra")
    mkpath(name)
    cd(name)
    l = 1
    #Recording data at initial timestep 
    writemacro(nnx,nny,ui,uj,rho,temp,printui,printuj,printrho,printtemp,solidnode,solidnodea,l)

    #Recording simulation parameters for reference if needed
    aparam = open("a_julia_m_parameters_$name.txt", "w") 
    write(aparam, "RE= $Re RA= $Ra PR= $Pr tcold= $tcold thot= $thot 
    \nnu= $nu mu= $mu alpha= $alphah beta= $betah u0= $Uo tauf= $tauf taut= $taut To= $To Fo= $Fo 
    \ndx= $dx dt= $dt nnx= $nnx nny= $nny Lo= $Lo ")
    close(aparam)
end 

#Begining of Main Progam
@timeit to "Main Loop" begin
    for n= 1 : nsteps

        @timeit to "stream" begin
            stream(nnx,nny,nnode,ndir,vpi,solidnode,solidstr,f1,f2,t1,t2,temp,thot,tcold,phi,phip,phim,invphip,omegat,uwxr,uwxl,uwxb,uwxt,uwyr,uwyl,uwyb,uwyt,cs2,dt,fbx,fby)
        end

        @timeit to "macrovar" begin
            macrovar(nnx,nny,nnode,f2,t2,rho,pres,cs2,fri,frj,ui,uj,temp,tmid,thot,tcold,betah,Fo,solidnode,uwxr,uwxl,uwxb,uwxt,uwyr,uwyl,uwyb,uwyt,dt)
        end

        @timeit to "vout" begin
            vout(nnx,nny,solidstr,rho,rhoold,ui,uiold,uj,ujold,temp,tempold,n)
        end

        @timeit to "writemacro" begin
            #Changes how often Macrovariable data is printed; default every 50000 steps
            if rem(n,500)==0
                writemacro(nnx,nny,ui,uj,rho,temp,printui,printuj,printrho,printtemp,solidnode,solidnodea,n)
            end
        end

        @timeit to "collide" begin
            collide_CM_thermal(tauf,taut,nnx,nny,nnode,solidnode,ui,uj,rho,f2,t2,cs2,temp,pres,fri,frj,tmid,Fo,betah)
        end

        @timeit to "switch" begin
            global fcopy = f1
            global f1 = f2 
            global f2 = fcopy

            global tcopy = t1
            global t1 = t2 
            global t2 = tcopy
        end

        #Error checks that will close simulation if 3 primary Macrovariable feilds reach specified convergence
        if errorx < 1e-16  && n > 50 
            if xconv == 1
                @timeit to "writemacro" begin
                    writemacro(nnx,nny,ui,uj,rho,temp,printui,printuj,printrho,printtemp,solidnode,solidnodea,n) 
                    errx = open("julia_m_errorx_$n.txt", "w") 
                    write(errx, "g2th = $errort, g2ui = $errorx, g2uj = $errory \n")
                    close(errx)
                    global xconv = 0
                end
            end
        end
        if errory < 1e-16  && n > 50
            if yconv == 1
                @timeit to "writemacro" begin
                    writemacro(nnx,nny,ui,uj,rho,temp,printui,printuj,printrho,printtemp,solidnode,solidnodea,n)
                    erry = open("julia_m_errory_$n.txt", "w") 
                    write(erry, "g2th = $errort, g2ui = $errorx, g2uj = $errory \n")
                    close(erry)
                    global yconv = 0
                end
            end
        end
        if errort < 1e-16 && n > 50
            if tconv == 1
                @timeit to "writemacro" begin
                    writemacro(nnx,nny,ui,uj,rho,temp,printui,printuj,printrho,printtemp,solidnode,solidnodea,n)
                    errt = open("julia_m_errort_$n.txt", "w") 
                    write(errt, "g2th = $errort, g2ui = $errorx, g2uj = $errory \n")
                    close(errt)
                    global tconv = 0
                end
            end
        end
        if (errort < 1e-16) && (errorx < 1e-16) && (errory < 1e-16)
            @timeit to "writemacro" begin
                writemacro(nnx,nny,ui,uj,rho,temp,printui,printuj,printrho,printtemp,solidnode,solidnodea,n)
                println("nstep ", n ,"|| g2th ", errort ,"|| g2ui ",errorx ,"|| g2uj ",errory)
                errf = open("julia_m_final_$n.txt", "w") 
                    write(errf, "g2th = $errort, g2ui = $errorx, g2uj = $errory \n")
                close(errf)
            break
            end
        end
    end 
end

#Displaying and saving Timing outputs
show(to)
open("timer.json","w") do io 
    JSON3.pretty(io,TimerOutputs.todict(to))
end

#Displaying results and directory location of simulation completed. 
println("")
println("")
println(name)
println("")
println("RE= ",Re," RA= ",Ra," PR= ", Pr )
println("")
println("dx= ",dx, " dt= ",dt," nu= ",nu, "mu= ",mu, "beta= ",betah," u0= ",Uo," tauf= ",tauf, "taut= ",taut," To= ",To," nnx= ",nnx," nny= ",nny)
println("")
println("")
println(@__DIR__)
println("")
println("")
cd()
