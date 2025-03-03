function fitCCEntEntScal(
    SvN_list::Vector{Float64};
    err::Vector{Float64}=0.0SvN_list,
    mincut::Int=1,
    pbc::Bool=false)

    # log of chord length / 6 for open boundary
    logChord(l, L) = @. log(sin(π * l /L))/6
    
    L = length(SvN_list) + 1

    # fit scaling
    lm(x,p) = @. p[1] * x + p[2]
    xdata = logChord([1:L-1;],L); #log.(sin.(π .* [1:L-1;] ./L))./6
    fit = curve_fit(lm, xdata[mincut:L-mincut], SvN_list[mincut:L-mincut], [0.5, 0.0])
    fitparam = fit.param
    cent = fitparam[1]
    cent_err = stderror(fit)[1]
    if pbc
        cent /= 2.0
        cent_err/= 2.0
    end
    @show cent ± cent_err

    # plot scaling
    fig = scatter(1:L-1, SvN_list, ylabel=L"S_{vN}", xlabel=L"l", frame=:box, yerror=err, label=false, lw=2, marker=:circle, xlims=(-1, L+1))
    plot!(1:L-1, fitparam[1] .* logChord([1:L-1;], L) .+ fitparam[2], label=false)

    # plot rescaled
    plot!(subplot=2, framestyle=:box,
    inset = (1, bbox(0.3, 0.2, 0.4, 0.45, :bottom)))
    xdata = LinRange(xdata[1], 0, 25)
    ydata = fit.param[1] * xdata .+ fit.param[2]

    if pbc
        scatter!(subplot=2, lw=2,
        log.(sin.(π .*[1:L-1;]./L)) ./3, SvN_list,
        xlabel=L"\frac{1}{3}\ln\sin(π l/L)",
        yerror=err, marker=:circle, label=false)
        plot!(subplot=2, 2xdata, ydata, lw=2,label=(@sprintf "c = %.2f ± %.2f" cent cent_err))
    else
        scatter!(subplot=2, lw=2,
        log.(sin.(π .*[1:L-1;]./L)) ./6, SvN_list,
        xlabel=L"\frac{1}{6}\ln\sin(π l/L)",
        yerror=err, marker=:circle, label=false)
        plot!(subplot=2, xdata, ydata, lw=2,label=(@sprintf "c = %.2f ± %.2f" cent cent_err))
    end
    return cent, fig
end

function new_fitCCEntEntScal(
    SvN_list::Vector{Float64};
    err::Vector{Float64}=0.0SvN_list,
    mincut::Int=1,
    pbc::Bool=false)

    # log of chord length / 6 for open boundary
    logChord(l, L) = @. log(sin(π * l /L))/6

    L = length(SvN_list) + 1
    mincut=mincut-1
    # fit scaling
    lm(x,p) = @. p[1] * x + p[2]
    xdata = logChord([1:L-1;],L); #log.(sin.(π .* [1:L-1;] ./L))./6
    fit = curve_fit(lm, xdata[setdiff(1:L-1, mincut:L-mincut)], SvN_list[setdiff(1:L-1, mincut:L-mincut)], [0.5, 0.0])
    fitparam = fit.param
    cent = fitparam[1]
    cent_err = stderror(fit)[1]
    if pbc
        cent /= 2.0
        cent_err/= 2.0
    end
    @show cent ± cent_err

    # plot rescaled
    plot!(subplot=2, framestyle=:box,
    inset = (1, bbox(0.3, 0.2, 0.4, 0.45, :bottom)))
    xdata = LinRange(xdata[1], 0, 25)
    ydata = fit.param[1] * xdata .+ fit.param[2]
    
    if pbc
        scatter!(subplot=2, lw=2,
        log.(sin.(π .*[1:L-1;]./L)) ./3, SvN_list,
        xlabel=L"\frac{1}{3}\ln\sin(π l/L)",
        yerror=err, marker=:circle, label=false)
        plot!(subplot=2, 2xdata, ydata, lw=2,label=(@sprintf "c = %.2f ± %.2f" cent cent_err))
    else
        scatter!(subplot=2, lw=2,
        log.(sin.(π .*[1:L-1;]./L)) ./6, SvN_list,
        xlabel=L"\frac{1}{6}\ln\sin(π l/L)",
        yerror=err, marker=:circle, label=false)
        plot!(subplot=2, xdata, ydata, lw=2,label=(@sprintf "c = %.2f ± %.2f" cent cent_err))
    end
    return cent, fig
end

function fitpage_curve(SvN_list::Vector{Float64}; err::Vector{Float64}=0.0SvN_list, mincut::Int=1)
    # Page curve formula is S(m,n)=lnm- m/2n, where m,n are the dimensions of the Hilbert space, and m<=n. If we take m=min(2^l,2^(L-l)), n=2^L, then S(l) = log(2)*min(l,L-l)- 2^(min(l,L-l))/(2^(L-min(l,L-2))*2)
    
    pagecurve(l, L) = @. log(2)*min(l,L-l) - 2^(min(l,L-l))/(2^(L-min(l,L-l))*2)
    L = length(SvN_list) + 1

    # fit scaling
    lm(x,p) = @. p[1] * x + p[2]
    xdata = pagecurve([1:L-1;],L); #log(2)*min(l,L-l)- 2^(min(l,L-l))/(2^(L-min(l,L-2))*2)
    fit = curve_fit(lm, xdata[mincut:L-mincut], SvN_list[mincut:L-mincut], [0.5, 0.0])
    fitparam = fit.param
    cent = fitparam[1]
    cent_err = stderror(fit)[1]
    @show cent ± cent_err

    # plot scaling
    fig = scatter(1:L-1, SvN_list, ylabel=L"S_{vN}", xlabel=L"l", frame=:box,
        yerror=err, label=false, lw=2, marker=:circle, xlims=(-1, L+1))
    plot!(1:L-1, fitparam[1] .* pagecurve([1:L-1;], L) .+ fitparam[2], label=false)
    
    # plot rescaled
    plot!(subplot=2, framestyle=:box,
        inset = (1, bbox(0.35, 0.15, 0.3, 0.35, :bottom)))

    scatter!(subplot=2, lw=2, xdata, SvN_list,
        xlabel=L"c(\ln2) l- 2^{2l-L-1}", yerror=err, marker=:circle, label=false)
    plot!(subplot=2, xdata, fit.param[1] * xdata .+ fit.param[2], lw=2,label=(@sprintf "c = %.2f ± %.2f" cent cent_err))

    return cent, fig

end

function fit_both(
    SvN_list::Vector{Float64},
    SvN_list2::Vector{Float64},
    err=0.0,
    mincut::Int=1,
    pbc::Bool=true)

    function fit_cc(SvN_list::Vector{Float64},mincut::Int=1,
        pbc::Bool=true)
        logChord(l, L) = @. log(sin(π * l /L))/6
       
        L = length(SvN_list) + 1

        # fit scaling
        lm(x,p) = @. p[1] * x + p[2]
        xdata = logChord([1:L-1;],L); 
        fit = curve_fit(lm, xdata[mincut:L-mincut], SvN_list[mincut:L-mincut], [0.5, 0.0])
        fitparam = fit.param
        cent = fitparam[1]
        cent_err = stderror(fit)[1]
        if pbc
            cent /= 2.0
            cent_err/= 2.0
        end
        @show cent ± cent_err

        return xdata,cent,cent_err, fitparam, fit
    end
    
    function fit_pagecurve(SvN_list::Vector{Float64},mincut::Int=1,
        pbc::Bool=true)
        lm(x,p) = @. p[1] * x + p[2]
        xdata = pagecurve([1:L-1;],L); #log(2)*min(l,L-l)- 2^(min(l,L-l))/(2^(L-min(l,L-2))*2)
        fit = curve_fit(lm, xdata[mincut:L-mincut], SvN_list[mincut:L-mincut], [0.5, 0.0])
        fitparam = fit.param
        cent = fitparam[1]
        cent_err = stderror(fit)[1]
        @show cent ± cent_err

        return xdata,cent,cent_err, fitparam, fit
    end

    xdata,cent,cent_err ,fitparam, fit=fit_cc(SvN_list, mincut, pbc)
    xdata2,cent2, cent_err2,fitparam2, fit2=fit_pagecurve(SvN_list2, mincut, pbc)
    L = length(SvN_list) + 1
    fig=scatter(1:L-1, SvN_list, ylabel=L"S_{vN}", xlabel=L"l", frame=:box, yerror=err, label=false, lw=2, marker=:circle, xlims=(-1, L+1),color="red",markersize=6)

    plot!(1:L-1, fitparam[1] .* logChord([1:L-1;], L) .+ fitparam[2], label=L"\frac{c}{3}\ln\sin(π l/L),  c=1.89 ",color="red", linestyle=:dash)

    scatter!(1:L-1, SvN_list2, ylabel=L"S_{vN}", xlabel=L"l", frame=:box, yerror=err, label=false, lw=2, marker=:circle, xlims=(-1, L+1),color=RGB(12/255, 159/255, 250/255),markersize=6)
    plot!( 1:L-1, fitparam2[1] .*  linearGrowthDecay([1:L-1;],L) .+ fitparam2[2], label=L"\alpha l-2^{2l-L-1}, \alpha=0.42", linestyle=:dash,color=RGB(12/255, 159/255, 250/255))

        # plot rescaled
        # plot!(subplot=2, framestyle=:box,
        # inset = (1, bbox(0.3, 0.2, 0.4, 0.45, :bottom)))


        # xdata = LinRange(xdata[1], 0, 25)
        # ydata = fit.param[1] * xdata .+ fit.param[2]
        # if pbc
        #     scatter!(subplot=2, lw=2,
        #     log.(sin.(π .*[1:L-1;]./L)) ./3, SvN_list,
        #     xlabel=L"\frac{1}{3}\ln\sin(π l/L)",
        #     yerror=err, marker=:circle, label=false)

        #     plot!(subplot=2, 2xdata, ydata, lw=2,label=(@sprintf "c = %.2f ± %.2f" cent cent_err))
        # else
        #     scatter!(subplot=2, lw=2,
        #     log.(sin.(π .*[1:L-1;]./L)) ./6, SvN_list,
        #     xlabel=L"\frac{1}{6}\ln\sin(π l/L)",
        #     yerror=err, marker=:circle, label=false)

        #     plot!(subplot=2, xdata, ydata, lw=2,label=(@sprintf "c = %.2f ± %.2f" cent cent_err))
        # end

    return cent, cent2, fig
end
