using LaTeXStrings, Printf, LsqFit, Measurements, Plots
 
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

function fitLpluslnL(SvN_list::Vector{Float64}; err::Vector{Float64}=0.0SvN_list, mincut::Int=1)
    # Now suppose such system has L + lnL SvN, like the majorana spin liquid.
    
    pagecurve(l, L) = @. log(2)*min(l,L-l) - 2^(min(l,L-l))/(2^(L-min(l,L-l))*2)
    logChord(l, L) = @. log(sin(π * l /L))/6
    llis= collect(1:length(SvN_list))
    L = length(SvN_list) + 1

    # fit scaling
    lm(x, p) = @. p[1] * pagecurve(x,L) + p[2] * logChord(x,L) + p[3]
    fit = curve_fit(lm, llis[mincut:L-mincut], SvN_list[mincut:L-mincut], [0.5, 0.5,0.0])
    params = coef(fit)
    error = stderror(fit)
    @show params, error

    # plot scaling
    fig = scatter(llis, SvN_list, ylabel=L"S_{vN}", xlabel=L"l", frame=:box,
        yerror=err, label=false, lw=2, marker=:circle, xlims=(-1, L+1))
    plot!(llis, lm(llis,params), label=false)
    
    # plot rescaled
    plot!(subplot=2, framestyle=:box,
        inset = (1, bbox(0.35, 0.15, 0.3, 0.35, :bottom)))

    xdata=pagecurve.(llis,L) + params[2]/params[1] .* logChord.(llis,L)
    scatter!(subplot=2, lw=2, xdata, SvN_list,
        xlabel=L"a(\ln2) l- 2^{2l-L-1} + b\frac{1}{3}\ln\sin(π l/L)", yerror=err, marker=:circle, label=false)
    plot!(subplot=2, xdata, params[1] * xdata .+ params[3], lw=2,label=(@sprintf "a=%.2f, b=%.2f" params[1] params[2]), 
    legend_background_color=nothing,
    legend_foreground_color=nothing,
    legendsize=1)

    return params[1:2], fig

end

function fit_both(
    llis::Vector{Int64},
    SvN_list::Vector{Float64},
    SvN_list2::Vector{Float64},
    err=0.0,
    mincut::Int=1,
    pbc::Bool=true)

    @assert length(llis) == length(SvN_list) "length of llis should be $(length(SvN_list)), but got $(length(llis))"
    
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
        L = length(SvN_list) + 1
        pagecurve(l, L) = @. log(2)*min(l,L-l) - 2^(min(l,L-l))/(2^(L-min(l,L-l))*2)
        xdata = pagecurve([1:L-1;],L); 
        lm(x,p) = @. p[1] * x + p[2]
       
        fit = curve_fit(lm, xdata[mincut:L-mincut], SvN_list[mincut:L-mincut], [0.5, 0.0])
        fitparam = fit.param
        cent = fitparam[1]
        cent_err = stderror(fit)[1]
        @show cent ± cent_err

        return xdata,cent,cent_err, fitparam, fit
    end

    logChord(l, L) = @. log(sin(π * l /L))/6
    pagecurve(l, L) = @. log(2)*min(l,L-l) - 2^(min(l,L-l))/(2^(L-min(l,L-l))*2)
    xdata,cent,cent_err ,fitparam, fit=fit_cc(SvN_list, mincut, pbc)
    xdata2,cent2, cent_err2,fitparam2, fit2=fit_pagecurve(SvN_list2, mincut, pbc)
    L = length(SvN_list) + 1

    fig=scatter(llis, SvN_list, ylabel=L"S_{vN}", xlabel=L"l", frame=:box, yerror=err, label=false, lw=2, marker=:circle, xlims=(-1, L+1),color="red",markersize=6, xticks=llis, legend_background_color=nothing,
    legend_foreground_color=nothing,
    )

    plot!(llis, fitparam[1] .* logChord(llis, L) .+ fitparam[2], label=latexstring("\\frac{c}{3}\\ln\\sin(π l/L),  c=$(round(cent,digits=2))"),color="red", linestyle=:dash)

    scatter!(llis, SvN_list2, ylabel=L"S_{vN}", xlabel=L"l", frame=:box, yerror=err, label=false, lw=2, marker=:circle, xlims=(-1, L+1),color=RGB(12/255, 159/255, 250/255),markersize=6)
    plot!(llis, fitparam2[1] .*  pagecurve([1:L-1;],L) .+ fitparam2[2], label=latexstring("\\alpha l-2^{2l-L-1}, \\alpha=$(round(fitparam2[1],digits=2))"), linestyle=:dash,color=RGB(12/255, 159/255, 250/255))

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

function ee_PXP_scaling_fig(N::Int64, state::Vector{ET},fit::String) where {ET}
    splitlis=Vector(1:N-1)
    EElis=ee_PXP_state(N, splitlis, state)

    if fit=="CC" 
        cent, fig=fitCCEntEntScal(EElis; mincut=1, pbc=true)
    end

    if fit=="Page"
        cent, fig=fitpage_curve(EElis; mincut=1)
    end

    if fit=="L+lnL"
        cent, fig=fitLpluslnL(EElis; mincut=1)
    end
    return cent, fig
end