# To run this from a terminal:
# > julia grok.jl

# To run this from within Julia:
# julia> include("grok.jl")


using Korg

# We should loop this but my Julia-foo is bad.
# Sun
photosphere_path = "/uufs/chpc.utah.edu/common/home/u6020307/grok/data/photospheres/marcs_krz/sun.krz"
metallicity = 0

function synthesize(atmosphere_path, linelist_path, lambdas, metallicity)
    println("Running from ", atmosphere_path, " with ", lambdas)
    atm = Korg.read_model_atmosphere(atmosphere_path)
    linelist = Korg.read_linelist(linelist_path)
    println("Synthesizing..")
    @time spectrum = Korg.synthesize(atm, linelist, lambdas; metallicity=metallicity)
    println("Done")
    return spectrum
end

println("Let's get to business")
@time spectrum = synthesize(photosphere_path, "/uufs/chpc.utah.edu/common/home/u6020307/grok/data/transitions/all-3930-3950.vald", Korg.air_to_vacuum(3930):0.01:Korg.air_to_vacuum(3950), metallicity)

# Save to disk.
open("executions/Sun_3930_3950_korg", "w") do fp
    for flux in spectrum.flux
        println(fp, flux)
    end
end

function synthesize_2(atmosphere_path, linelist_path_a, linelist_path_b, lambdas, metallicity)
    println("Running from ", atmosphere_path, " with ", lambdas)
    atm = Korg.read_model_atmosphere(atmosphere_path)
    linelist_a = Korg.read_linelist(linelist_path_a)
    linelist_b = Korg.read_linelist(linelist_path_b)
    linelist = append!(linelist_a, linelist_b)
    println("Synthesizing..")
    @time spectrum = Korg.synthesize(atm, linelist, lambdas; metallicity=metallicity)
    println("Done")
    return spectrum
end

@time spectrum = synthesize_2(
    photosphere_path, 
    "/uufs/chpc.utah.edu/common/home/u6020307/grok/data/transitions/all-5160-5176.vald",
    "/uufs/chpc.utah.edu/common/home/u6020307/grok/data/transitions/all-5176-5190.vald",
    Korg.air_to_vacuum(5160):0.01:Korg.air_to_vacuum(5190),
    metallicity
)
# Save to disk.
open("executions/Sun_5160_5190_korg", "w") do fp
    for flux in spectrum.flux
        println(fp, flux)
    end
end

@time spectrum = synthesize_2(
    photosphere_path, 
    "/uufs/chpc.utah.edu/common/home/u6020307/grok/data/transitions/all-6540-6559.vald", 
    "/uufs/chpc.utah.edu/common/home/u6020307/grok/data/transitions/all-6559-6578.vald", 
    Korg.air_to_vacuum(6540):0.01:Korg.air_to_vacuum(6578),
    metallicity
)
# Save to disk.
open("executions/Sun_6540_6578_korg", "w") do fp
    for flux in spectrum.flux
        println(fp, flux)
    end
end



# Arcturus
photosphere_path = "/uufs/chpc.utah.edu/common/home/u6020307/s4250_g+1.5_m1.0_t02_st_z-0.50_a+0.20_c+0.00_n+0.00_o+0.20_r+0.00_s+0.00.krz"
metallicity = -0.53

function synthesize(atmosphere_path, linelist_path, lambdas, metallicity)
    println("Running from ", atmosphere_path, " with ", lambdas)
    atm = Korg.read_model_atmosphere(atmosphere_path)
    linelist = Korg.read_linelist(linelist_path)
    println("Synthesizing..")
    @time spectrum = Korg.synthesize(atm, linelist, lambdas; metallicity=metallicity)
    println("Done")
    return spectrum
end

println("Doing initial allocation")
@time spectrum = synthesize(photosphere_path, "/uufs/chpc.utah.edu/common/home/u6020307/grok/data/transitions/all-3930-3950.vald", Korg.air_to_vacuum(3930):0.01:Korg.air_to_vacuum(3950), metallicity)

println("Let's get to business")
@time spectrum = synthesize(photosphere_path, "/uufs/chpc.utah.edu/common/home/u6020307/grok/data/transitions/all-3930-3950.vald", Korg.air_to_vacuum(3930):0.01:Korg.air_to_vacuum(3950), metallicity)

# Save to disk.
open("executions/Arcturus_3930_3950_korg", "w") do fp
    for flux in spectrum.flux
        println(fp, flux)
    end
end

function synthesize_2(atmosphere_path, linelist_path_a, linelist_path_b, lambdas, metallicity)
    println("Running from ", atmosphere_path, " with ", lambdas)
    atm = Korg.read_model_atmosphere(atmosphere_path)
    linelist_a = Korg.read_linelist(linelist_path_a)
    linelist_b = Korg.read_linelist(linelist_path_b)
    linelist = append!(linelist_a, linelist_b)
    println("Synthesizing..")
    @time spectrum = Korg.synthesize(atm, linelist, lambdas; metallicity=metallicity)
    println("Done")
    return spectrum
end

@time spectrum = synthesize_2(
    photosphere_path, 
    "/uufs/chpc.utah.edu/common/home/u6020307/grok/data/transitions/all-5160-5176.vald",
    "/uufs/chpc.utah.edu/common/home/u6020307/grok/data/transitions/all-5176-5190.vald",
    Korg.air_to_vacuum(5160):0.01:Korg.air_to_vacuum(5190),
    metallicity
)
# Save to disk.
open("executions/Arcturus_5160_5190_korg", "w") do fp
    for flux in spectrum.flux
        println(fp, flux)
    end
end

@time spectrum = synthesize_2(
    photosphere_path, 
    "/uufs/chpc.utah.edu/common/home/u6020307/grok/data/transitions/all-6540-6559.vald", 
    "/uufs/chpc.utah.edu/common/home/u6020307/grok/data/transitions/all-6559-6578.vald", 
    Korg.air_to_vacuum(6540):0.01:Korg.air_to_vacuum(6578),
    metallicity
)
# Save to disk.
open("executions/Arcturus_6540_6578_korg", "w") do fp
    for flux in spectrum.flux
        println(fp, flux)
    end
end




#@time begin
#    print("Running 3660 - 3680 for STAR")
#    atm = Korg.read_model_atmosphere(photosphere_path)
#    linelist = Korg.read_linelist("/uufs/chpc.utah.edu/common/home/u6020307/grok/data/transitions/all-3660-3680.vald")
#    synthesis = Korg.synthesize(atm, linelist, 3660:0.01:3680)
#end

# ERROR: LoadError: DomainError with 3649.0:
#The wavelength must lie in the interval [3847 Å, 25000 Å]

