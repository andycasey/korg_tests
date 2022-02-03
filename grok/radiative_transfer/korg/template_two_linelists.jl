using Pkg

Pkg.status("Korg")

using Korg

function synthesize(atmosphere_path, linelist_path_0, linelist_path_1, metallicity)
    println("Running with ", atmosphere_path, " and ", linelist_path_0, " and ", linelist_path_1)
    atm = Korg.read_model_atmosphere(atmosphere_path)
    linelist_1 = Korg.read_linelist(linelist_path_0, format="{korg_read_transitions_format}")
    linelist_2 = Korg.read_linelist(linelist_path_1, format="{korg_read_transitions_format}")
    linelist = append!(linelist_1, linelist_2)
    println("Synthesizing..")
    @time spectrum = Korg.synthesize(atm, linelist, {lambda_vacuum_min:.2f}, {lambda_vacuum_max:.2f}; metallicity=metallicity, hydrogen_lines={hydrogen_lines}, vmic={microturbulence:.2f})
    println("Done")
    return spectrum
end

function synthesize_one(atmosphere_path, linelist_path, metallicity)
    println("Running with ", atmosphere_path, " and ", linelist_path)
    atm = Korg.read_model_atmosphere(atmosphere_path)
    linelist = Korg.read_linelist(linelist_path, format="{korg_read_transitions_format}")
    println("Synthesizing..")
    @time spectrum = Korg.synthesize(atm, linelist, {lambda_vacuum_min:.2f}, {lambda_vacuum_max:.2f}; metallicity=metallicity, hydrogen_lines={hydrogen_lines}, vmic={microturbulence:.2f})
    println("Done")
    return spectrum
end

println("Going once.")
@time spectrum = synthesize("{atmosphere_path}", "{linelist_path_0}", "{linelist_path_1}", {fake_metallicity:.2f})

println("Going twice..")
@time spectrum = synthesize("{atmosphere_path}", "{linelist_path_0}", "{linelist_path_1}", {metallicity:.2f})

# Save to disk.
println("Going three times...")
open("spectrum.out", "w") do fp
    for (wl, flux) in zip(spectrum.wavelengths, spectrum.flux)
        println(fp, wl, " ", flux)
    end
end
println("Sold!")

# Now do continuum.
atm = Korg.read_model_atmosphere("{atmosphere_path}")
continuum = Korg.synthesize(atm, [], {lambda_vacuum_min:.2f}, {lambda_vacuum_max:.2f}; metallicity={metallicity:.2f}, hydrogen_lines=false, vmic={microturbulence:.2f})
open("continuum.out", "w") do fp
    for (wl, flux) in zip(continuum.wavelengths, continuum.flux)
        println(fp, wl, " ", flux)
    end
end