# =============================
#     Program: plot_lib.jl
#     Aim: functions for plotting
#     Author: Sergio Pérez-Montero, 2022.11.11
# =============================

function plot_spectrum(x, d::Any, f::Any, G::Any, vrs::Any, clrmp, plotpath::String; fntsz=nothing, fancy=false)
    ## First determine plot parameters
    # -- number of rows and columns
    nrows, ncols = length(vrs), 2

    # -- figure size
    fgsz = (1500*ncols, 500*nrows)

    # -- check fontsize
    isnothing(fntsz) && (fntsz = 0.01 * sqrt(fgsz[1]^2 + fgsz[2]^2))
    fontsize_theme = Theme(font="Dejavu Serif", fontsize=fntsz)
    set_theme!(fontsize_theme)

    # Now, define the figure and static elements
    fig = Figure(resolution=fgsz)

    ## Plotting
    # -- Plot variables
    for i in 1:nrows
        if i == 1
            ax = Axis(fig[i, 1], title = "AMOD variables", xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Time (kyr)", ylabel=vrs[i]*" ( "*d[i].attrib["units"]*")", xgridcolor=:darkgrey, ygridcolor=:darkgrey)
            ax_spect = Axis(fig[i, 2], title="Normalized Power density", xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Period"*" ("*x.attrib["units"]*")", xgridcolor=:darkgrey, ygridcolor=:darkgrey)
        else
            ax = Axis(fig[i, 1], titlesize=0.7 * fntsz, xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Time (kyr)", ylabel=vrs[i]*" ("*d[i].attrib["units"]*")", xgridcolor=:darkgrey, ygridcolor=:darkgrey)
            ax_spect = Axis(fig[i, 2], xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Period"*" (k"*x.attrib["units"]*")", xgridcolor=:darkgrey, ygridcolor=:darkgrey)
        end
        if i < nrows
            hidexdecorations!(ax, grid=false)
            hidexdecorations!(ax_spect, grid=false)
        end
        
        update_theme!()

        maxdi, mindi = maximum(d[i][:]), minimum(d[i][:])
        maxi = max(abs(maxdi), abs(mindi))
        if var(d[i]) < 1e-1
            lines!(ax, x, d[i], linewidth=3, color = x, colormap = clrmp[i])
        else
            lines!(ax, x, d[i], linewidth=3, color = d[i], colormap = clrmp[i], colorrange = (-maxi, maxi))
        end

        lines!(ax_spect, f[i], G[i], color=:black, linewidth=3)
        if fancy
            band!(ax_spect, f[i], 0.0, G[i], linewidth=2, color = G[i], colormap = :berlin)
        else
            band!(ax_spect, f[i], 0.0, G[i], linewidth=2, color = :darkred)
        end

        xlims!(ax, (x[1], x[end]))
        xlims!(ax_spect, (1/500e3, 1/21e3))

        xlen = length(x)
        if mod(xlen, 2) == 0
            xstep = Int(xlen / 10)
        else
            xstep = Int((xlen - 1) / 10)
        end 
        ax.xticks = x[1:xstep:end]
        ax.xtickformat = k -> string.(k/1000)

        ax_spect.xticks = [1/100e3, 1/41e3, 1/23e3]
        ax_spect.xtickformat = k -> string.(Int.(ceil.(1 ./ k)/1000)) 
    end

    # Resizing
    colsize!(fig.layout, 2, Relative(0.33))
    resize_to_layout!(fig)

    # Saving
    save(plotpath, fig)
end

function plot_amod(experiment="test_default", vars = ["ins_norm", "SMB", "H", "Hsed"])
    ## Load output data 
    data, time = load_nc(amod_path * "/output/" * experiment * "/amod.nc", vars)

    ## Calculate spectra
    #window = 50
    G_data, freqs_data = [], []
    for v in 1:length(vars)
        new_data = copy(data[v])

        # -- spectrum
        if false
            # ---- fft
            F, N = fft(new_data), length(new_data)
            Nmax = Int(ceil(N/2))
            kvec = 1:Nmax
            G = abs.(F)*2/N
            G = sqrt.(G[1:Nmax])
            freq = kvec ./ (N*(time[2] - time[1])) 
            G, freq = G[freq .>= 1/150.5e3], freq[freq .>= 1/150e3] # we eliminate values above and below Milankovitch cycles
            G, freq = G[freq .<= 1/22e3], freq[freq .<= 1/22e3]
            G = G ./ sum(G)
            if var(new_data) < 1e-1
                G = zeros(length(G))
            end
        else
            # ---- blackman tuckey
            N = length(new_data)
            Nmax = Int(ceil(N/2))
            P = periodogram(new_data, onesided=false, fs=1/(time[2] - time[1]), window=blackman(N))
            G, freq = P.power, P.freq
            G, freq = G[1:Nmax] .* 2, freq[1:Nmax]
            G = sqrt.(G) / (N/2)
            G, freq = G[freq .>= 1/150.5e3], freq[freq .>= 1/150e3] # we eliminate values above and below Milankovitch cycles
            G, freq = G[freq .<= 1/22e3], freq[freq .<= 1/22e3]
            G = G ./ sum(G)
            if var(new_data) < 1e-1
                G = zeros(length(G))
            end
        end
        push!(G_data, G)
        push!(freqs_data, freq)  # save frequencies

        # -- time vs freq map WORK IN PROGRESS -- spm 2022.11.18
        # calc_TimeFreqArr(new_data, window)

    end

    ## Plot
    clrmp_dict = Dict("time"=>:grays1,
        "H" => cgrad([:royalblue, :royalblue4]),#cgrad([:snow4, :royalblue, :royalblue4, :navy]),
        "Hsed" => cgrad([:firebrick, :maroon]),#cgrad([:maroon, :black]),
        "T" => :lajolla,
        "A" => :lajolla,
        "T_sl" => cgrad([:purple, :purple]),
        "TMB" => cgrad(:seismic, [0.0, 0.45, 0.55, 1.0], rev = true),
        "SMB" => cgrad(:seismic, [0.0, 0.45, 0.55, 1.0], rev = true),
        "Z" => :dense,
        "B" => :dense,
        "M" => :amp,
        "Acc" => cgrad(:ice, rev = true),
        "U_d" => :corkO,
        "U_b" => :corkO,
        "U" => :corkO,
        "T_surf" => :lajolla,
        "tau_b" => :thermal,
        "tau_d" => :thermal,
        "Q_dif" => :lajolla,
        "Q_difup" => :lajolla,
        "Q_difdown" => :lajolla,
        "Q_drag" => :lajolla,
        "alpha" => :lajolla,
        "Q_adv" => :lajolla,
        "fstream" => :grays1,
        "fstream_ref" => :grays1,
        "Hdot" => :ice,
        "Hseddot" => :copper,
        "Bdot" => :dense,
        "Tdot" => :lajolla,
        "fstreamdot" => :grays1,
        "ins" => cgrad([:black, :black]),#cgrad(:seismic, [0.0, 0.49, 0.51, 1.0]),
        "ins_norm" => cgrad([:black, :black]),#cgrad(:seismic, [0.0, 0.49, 0.51, 1.0]),
        "co2" => :thermal,
        "P" => :dense,
        "exc" => :thermal,
        "long_peri" => :thermal,
        "obl" => :thermal
    )

    colors = [clrmp_dict[i] for i in vars]
    plot_spectrum(time, data, freqs_data, G_data, vars, colors, amod_path * "/output/" * experiment * "/" * "amod_results.png");

end