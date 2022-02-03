using Pkg

Pkg.status("Korg")

using Korg

function synthesize(atmosphere_path, linelist_path, metallicity)
    println("Running with ", atmosphere_path, " and ", linelist_path)
    atm = Korg.read_model_atmosphere(atmosphere_path)
    linelist = Korg.read_linelist(linelist_path, format="{korg_read_transitions_format}")
    println("Synthesizing..")
    @time spectrum = Korg.synthesize(atm, linelist, {lambda_vacuum_min}, {lambda_vacuum_max}; metallicity=metallicity, hydrogen_lines={hydrogen_lines}, vmic={microturbulence:.2f})
    println("Done")
    return spectrum
end

println("Going once.")
@time spectrum = synthesize("{atmosphere_path}", "{linelist_path}", {fake_metallicity:.2f})

println("Going twice..")
@time spectrum = synthesize("{atmosphere_path}", "{linelist_path}", {metallicity:.2f})


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
continuum = Korg.synthesize(atm, [], {lambda_vacuum_min}, {lambda_vacuum_max}; metallicity={metallicity:.2f}, hydrogen_lines=false, vmic={microturbulence:.2f})
open("continuum.out", "w") do fp
    for (wl, flux) in zip(continuum.wavelengths, continuum.flux)
        println(fp, wl, " ", flux)
    end
end