function calcAnalytical(prms::params{R,I}, rstate, βstate)  where {R<:Real, I<:Integer}

    @unpack  δ_th, α, θ, θ0, γ, W_bar, κ_r, κ_β, σ2_r, σ2_β, r_bar, τ_max, τ, θvec, θ0vec, αvec, y_idxs, Δτ, Δt_y = prms

    # find initial condition
    output1 = nlsolve(x -> calcM0(x, prms), [0.5], iterations=1000000)
    hbk = output1.zero[1]
    ZAB = σ2_r*(γ/W_bar)*θ*Hint(κ_r,δ_th[1],δ_th[2],τ_max)/(1 + σ2_r*(γ/W_bar)*α*(Fint(κ_r,κ_r,δ_th[1],τ_max) - Fint(κ_r,hbk[1],δ_th[1],τ_max))/(hbk[1] - κ_r))

    M0 = zeros(4)
    M0[1] = κ_r + σ2_r*(γ/W_bar)*α*Fint(κ_r,κ_r,δ_th[1],τ_max)
    M0[2] = σ2_β*(γ/W_bar)*α*ZAB*(Fint(κ_r,κ_r,δ_th[1],τ_max) - Fint(κ_r,hbk[1],δ_th[1],τ_max))/(hbk[1] - κ_r)
    M0[3] = -σ2_r*((γ/W_bar)*θ*Hint(κ_r,δ_th[1],δ_th[2],τ_max) - (γ/W_bar)*α*ZAB*(Fint(κ_r,κ_r,δ_th[1],τ_max) - Fint(κ_r,hbk[1],δ_th[1],τ_max))/(hbk[1] - κ_r))
    M0[4] = κ_β - σ2_β*ZAB*((γ/W_bar)*θ*(Hint(κ_r,δ_th[1],δ_th[2],τ_max) - Hint(hbk[1],δ_th[1],δ_th[2],τ_max)) - (γ/W_bar)*α*ZAB*(Fint(κ_r,κ_r,δ_th[1],τ_max) + Fint(hbk[1],hbk[1],δ_th[1],τ_max) - 2*Fint(κ_r,hbk[1],δ_th[1],τ_max))/(hbk[1] - κ_r))/(hbk[1] - κ_r)

    # find the M matrix
    output2 = nlsolve(x -> calcMvec(x, prms), M0, iterations=1000000)
    Mvec = output2.zero

    # reorganize as matrix
    Mmat = [Mvec[1] Mvec[2]; Mvec[3] Mvec[4]]
       
    # compute eigenvalues
    L = eigvals(Mmat)
    L2 = L[1]
    L1 = L[2]

    # compute eigenvectors
    ϕ1 = -(Mvec[1]-L1)/Mvec[3]
    ϕ2 = -(Mvec[1]-L2)/Mvec[3]

    # compute coefs in summation across factors
    c_r1 = -ϕ2/(ϕ1-ϕ2)
    c_r2 = ϕ1/(ϕ1-ϕ2)
    c_β1 = 1/(ϕ1-ϕ2)
    c_β2 = -1/(ϕ1-ϕ2)

    # compute A elements
    Ar = c_r1.*(1 .- exp.(-L1.*τ))./L1 + c_r2.*(1 .- exp.(-L2.*τ))./L2
    Aβ = c_β1.*(1 .- exp.(-L1.*τ))./L1 + c_β2.*(1 .- exp.(-L2.*τ))./L2
    
    # compute χ elements
    
    I_rth0 = c_r1*Hint(L1,δ_th[1],δ_th[2],τ_max) + c_r2*Hint(L2,δ_th[1],δ_th[2],τ_max)
    I_βth0 = c_β1*Hint(L1,δ_th[1],δ_th[2],τ_max) + c_β2*Hint(L2,δ_th[1],δ_th[2],τ_max)

    II_rr = c_r1^2*Lint(L1,L1,δ_th[1],τ_max) + c_r1*c_r2*(Lint(L1,L2,δ_th[1],τ_max)+Lint(L2,L1,δ_th[1],τ_max)) + c_r2^2*Lint(L2,L2,δ_th[1],τ_max)
    II_rβ = c_r1*c_β1*Lint(L1,L1,δ_th[1],τ_max) + c_r1*c_β2*Lint(L1,L2,δ_th[1],τ_max) + c_r2*c_β1*Lint(L2,L1,δ_th[1],τ_max) + c_r2*c_β2*Lint(L2,L2,δ_th[1],τ_max)
    II_βr = c_r1*c_β1*Lint(L1,L1,δ_th[1],τ_max) + c_r1*c_β2*Lint(L2,L1,δ_th[1],τ_max) + c_r2*c_β1*Lint(L1,L2,δ_th[1],τ_max) + c_r2*c_β2*Lint(L2,L2,δ_th[1],τ_max)
    II_ββ = c_β1^2*Lint(L1,L1,δ_th[1],τ_max) + c_β1*c_β2*(Lint(L1,L2,δ_th[1],τ_max) + Lint(L2,L1,δ_th[1],τ_max)) + c_β2^2*Lint(L2,L2,δ_th[1],τ_max)

    II_r2r = c_r1*(c_r1^2*Kint(L1,L1,L1,δ_th[1],τ_max)+c_r2^2*Kint(L2,L2,L1,δ_th[1],τ_max)+2*c_r1*c_r2*Kint(L1,L2,L1,δ_th[1],τ_max))+c_r2*(c_r1^2*Kint(L1,L1,L2,δ_th[1],τ_max)+c_r2^2*Kint(L2,L2,L2,δ_th[1],τ_max)+2*c_r1*c_r2*Kint(L1,L2,L2,δ_th[1],τ_max))
    II_r2β = c_β1*(c_r1^2*Kint(L1,L1,L1,δ_th[1],τ_max)+c_r2^2*Kint(L2,L2,L1,δ_th[1],τ_max)+2*c_r1*c_r2*Kint(L1,L2,L1,δ_th[1],τ_max))+c_β2*(c_r1^2*Kint(L1,L1,L2,δ_th[1],τ_max)+c_r2^2*Kint(L2,L2,L2,δ_th[1],τ_max)+2*c_r1*c_r2*Kint(L1,L2,L2,δ_th[1],τ_max))
    II_β2r = c_r1*(c_β1^2*Kint(L1,L1,L1,δ_th[1],τ_max)+c_β2^2*Kint(L2,L2,L1,δ_th[1],τ_max)+2*c_β1*c_β2*Kint(L1,L2,L1,δ_th[1],τ_max))+c_r2*(c_β1^2*Kint(L1,L1,L2,δ_th[1],τ_max)+c_β2^2*Kint(L2,L2,L2,δ_th[1],τ_max)+2*c_β1*c_β2*Kint(L1,L2,L2,δ_th[1],τ_max))
    II_β2β = c_β1*(c_β1^2*Kint(L1,L1,L1,δ_th[1],τ_max)+c_β2^2*Kint(L2,L2,L1,δ_th[1],τ_max)+2*c_β1*c_β2*Kint(L1,L2,L1,δ_th[1],τ_max))+c_β2*(c_β1^2*Kint(L1,L1,L2,δ_th[1],τ_max)+c_β2^2*Kint(L2,L2,L2,δ_th[1],τ_max)+2*c_β1*c_β2*Kint(L1,L2,L2,δ_th[1],τ_max))

    Convr = 0.5*(σ2_r*II_r2r + σ2_β*II_β2r)
    Convb = 0.5*(σ2_r*II_r2β + σ2_β*II_β2β)

    chir = ((κ_r*r_bar + σ2_r*(γ/W_bar)*(θ0*I_rth0 + α*Convr))*(1 + (γ/W_bar)*α*σ2_β*II_ββ) - σ2_β*(γ/W_bar)*(θ0*I_βth0 + α*Convb)*(γ/W_bar)*α*σ2_r*II_βr)/((1 + (γ/W_bar)*α*σ2_r*II_rr)*(1 + (γ/W_bar)*α*σ2_β*II_ββ) - ((γ/W_bar)*α)^2*σ2_r*σ2_β*II_rβ*II_βr)
    chib = (σ2_β*(γ/W_bar)*(θ0*I_βth0 + α*Convb)*(1 + (γ/W_bar)*α*σ2_r*II_rr)-(κ_r*r_bar + σ2_r*(γ/W_bar)*(θ0*I_rth0 + α*Convr))*(γ/W_bar)*α*σ2_β*II_rβ)/((1 + (γ/W_bar)*α*σ2_r*II_rr)*(1 + (γ/W_bar)*α*σ2_β*II_ββ) - ((γ/W_bar)*α)^2*σ2_r*σ2_β*II_rβ*II_βr)

    # compute C

    Cr = chir.*(c_r1.*(τ - (1 .- exp.(-L1.*τ))./L1)./L1 + c_r2.*(τ - (1 .- exp.(-L2.*τ))./L2)./L2)
    Cb = chib.*(c_β1.*(τ - (1 .- exp.(-L1.*τ))./L1)./L1 + c_β2.*(τ - (1 .- exp.(-L2.*τ))./L2)./L2)
    Crr = 0.5.*σ2_r.*((c_r1^2).*(τ .- 2*(1 .- exp.(-L1*τ))./L1 .+ (1 .- exp.(-2*L1*τ))./(2*L1))./L1^2 .+ c_r2^2*(τ .- 2*(1 .- exp.(-L2*τ))./L2 .+ (1 .- exp.(-2*L2*τ))./(2*L2))./L2^2 .+ 2*c_r1*c_r2*(τ .- (1 .- exp.(-L1*τ))./L1 .- (1 .- exp.(-L2*τ))./L2 .+ (1 .- exp.(-(L1+L2)*τ))./(L1+L2))./(L1*L2))
    Cbb = 0.5.*σ2_β.*((c_β1^2).*(τ .- 2*(1 .- exp.(-L1*τ))./L1 .+ (1 .- exp.(-2*L1*τ))./(2*L1))./L1^2 .+ c_β2^2*(τ .- 2*(1 .- exp.(-L2*τ))./L2 .+ (1 .- exp.(-2*L2*τ))./(2*L2))./L2^2 .+ 2*c_β1*c_β2*(τ .- (1 .- exp.(-L1*τ))./L1 .- (1 .- exp.(-L2*τ))./L2 .+ (1 .- exp.(-(L1+L2)*τ))./(L1+L2))./(L1*L2))

    C = Cr+Cb-Crr-Cbb

    # Truncate to exclude first entry
    Ar_tr = Ar[2:end]
    Aβ_tr = Aβ[2:end]
    C_tr  = C[2:end]
    τ_tr  = τ[2:end]

    # Compute moments
    y_bar = 100*(Ar_tr*r_bar + C_tr)./τ_tr
    σy  = 100*sqrt.((Ar_tr.^2)*σ2_r/(2*κ_r) + (Aβ_tr.^2)*σ2_β/(2*κ_β))./τ_tr
    σΔy = 100*((2*((σ2_r/(2*κ_r))*(1-exp(-κ_r))*(Ar_tr.^2)+(σ2_β/(2*κ_β))*(1-exp(-κ_β))*(Aβ_tr.^2))).^0.5)./τ_tr
    ρy  = (100^2)*((σ2_r/(2*κ_r))*(Ar_tr./τ_tr)*(Ar_tr./τ_tr)'+(σ2_β/(2*κ_β))*(Aβ_tr./τ_tr)*(Aβ_tr./τ_tr)')./(σy*σy')
    ρΔy = (100^2)*2*((σ2_r/(2*κ_r))*(1-exp(-κ_r))*(Ar_tr./τ_tr)*(Ar_tr./τ_tr)'+(σ2_β/(2*κ_β))*(1-exp(-κ_β))*(Aβ_tr./τ_tr)*(Aβ_tr./τ_tr)')./(σΔy*σΔy')

    # FB
    σ2_r_unc = σ2_r/(2*κ_r)
    σ2_β_unc = σ2_β/(2*κ_β)
    NFBr = (Ar_tr[y_idxs[2:30].-1] - Ar_tr[y_idxs[1:29].-1]*exp(-κ_r) .- Ar_tr[y_idxs[1]-1]).*(Ar_tr[y_idxs[2:30].-1] - Ar_tr[y_idxs[1:29].-1] .- Ar_tr[y_idxs[1]-1])
    NFBβ = (Aβ_tr[y_idxs[2:30].-1] - Aβ_tr[y_idxs[1:29].-1]*exp(-κ_β) .- Aβ_tr[y_idxs[1]-1]).*(Aβ_tr[y_idxs[2:30].-1] - Aβ_tr[y_idxs[1:29].-1] .- Aβ_tr[y_idxs[1]-1])
    FBn = NFBr*σ2_r_unc + NFBβ*σ2_β_unc
    FBd = ((Ar_tr[y_idxs[2:30].-1] - Ar_tr[y_idxs[1:29].-1] .- Ar_tr[y_idxs[1]-1]).^2)*σ2_r_unc + ((Aβ_tr[y_idxs[2:30].-1] - Aβ_tr[y_idxs[1:29].-1] .- Aβ_tr[y_idxs[1]-1]).^2)*σ2_β_unc
    FB = FBn./FBd

    amoments = zeros(10)
    amoments[1] = y_bar[y_idxs[1]-1]
    amoments[2] = y_bar[y_idxs[10]-1] - y_bar[y_idxs[1]-1]
    amoments[3] = y_bar[y_idxs[20]-1] - y_bar[y_idxs[1]-1]
    amoments[4] = FB[9]
    amoments[5] = σy[y_idxs[5]-1]
    amoments[6] = σy[y_idxs[10]-1]
    amoments[7] = σy[y_idxs[20]-1]
    amoments[8] = σΔy[y_idxs[20]-1]
    #amoments[6] = mean(σy[y_idxs[1:30].-1])
    #amoments[7] = mean(σΔy[y_idxs[1:30].-1])
    #amoments[8] = mean(ρΔy[y_idxs[1]-1,y_idxs[1:30].-1])
    
    #amoments[10] = FB[19]

    if rstate!=[]
        aterm = 100*(Ar[y_idxs].*rstate + Aβ[y_idxs].*βstate + C[y_idxs])./τ[y_idxs]
    else
        aterm = []
    end

    return amoments, aterm

end

function calcMvec(Mvec::Vector{R}, prms::params{R,I}) where {R<:Real, I<:Integer}


    @unpack  δ_th, α, γ, W_bar, θ, θ0, κ_r, κ_β, σ2_r, σ2_β, τ_max = prms

    # reorganize as matrix
    Mmat = [Mvec[1] Mvec[2]; Mvec[3] Mvec[4]]

    # compute eigenvalues
    L = eigvals(Mmat)
    L2 = L[1]
    L1 = L[2]

    # compute eigenvectors
    ϕ1 = -(Mvec[1]-L1)/Mvec[3]
    ϕ2 = -(Mvec[1]-L2)/Mvec[3]

    # compute coefs in summation across factors
    c_r1 = -ϕ2/(ϕ1-ϕ2)
    c_r2 = ϕ1/(ϕ1-ϕ2)
    c_β1 = 1/(ϕ1-ϕ2)
    c_β2 = -1/(ϕ1-ϕ2)

    I_r = c_r1^2*Fint(L1,L1,δ_th[1],τ_max) + c_r2^2*Fint(L2,L2,δ_th[1],τ_max) + 2*c_r1*c_r2*Fint(L1,L2,δ_th[1],τ_max)
    I_β = c_β1^2*Fint(L1,L1,δ_th[1],τ_max) + c_β2^2*Fint(L2,L2,δ_th[1],τ_max) + 2*c_β1*c_β2*Fint(L1,L2,δ_th[1],τ_max)
    I_rβ = c_r1*c_β1*Fint(L1,L1,δ_th[1],τ_max) + c_r2*c_β2*Fint(L2,L2,δ_th[1],τ_max) + (c_r1*c_β2+c_r2*c_β1)*Fint(L1,L2,δ_th[1],τ_max)
    I_rth = c_r1*Hint(L1,δ_th[1],δ_th[2],τ_max) + c_r2*Hint(L2,δ_th[1],δ_th[2],τ_max)
    I_βth = c_β1*Hint(L1,δ_th[1],δ_th[2],τ_max) + c_β2*Hint(L2,δ_th[1],δ_th[2],τ_max)

    # reconstruct M vector
    Mvec_new = zeros(4)
    Mvec_new[1] = κ_r + (γ/W_bar)*α*σ2_r*I_r
    Mvec_new[2] = (γ/W_bar)*α*σ2_β*I_rβ
    Mvec_new[3] = -σ2_r*(γ/W_bar)*(θ*I_rth - α*I_rβ)
    Mvec_new[4] = κ_β - σ2_β*(γ/W_bar)*(θ*I_βth - α*I_β)

    # output
    Mdiff = Mvec - Mvec_new

    return Mdiff

end

# Compute integrals

function Fint(l1, l2, δ, τ_max)
        
    F = 1/(l1*l2)*((1-exp(-δ*τ_max))/δ-(1-exp(-(l1+δ)*τ_max))/(l1+δ)-(1-exp(-(l2+δ)*τ_max))/(l2+δ)+(1-exp(-(l1+l2+δ)*τ_max))/(l1+l2+δ))

end

function Gint(l, δ, τ_max)
        
    G = 1/l*((1-exp(-δ*τ_max)-δ*τ_max*exp(-δ*τ_max))/δ^2-(1-exp(-(δ+l)*τ_max)-(δ+l)*τ_max*exp(-(δ+l)*τ_max))/(δ+l)^2)

end

function Hint(l, δ_al, δ_th, τ_max)
        
    H = 1/l*((1-exp(-δ_al*τ_max))/δ_al-(1-exp(-(δ_al+l)*τ_max))/(δ_al+l)-(1-exp(-δ_th*τ_max))/δ_th+(1-exp(-(δ_th+l)*τ_max))/(δ_th+l))

end

function Kint(l1, l2, l3, δ, τ_max)
        
    K = 1/(l1*l2)*(Gint(l3,δ,τ_max)-Fint(l1,l3,δ,τ_max)-Fint(l2,l3,δ,τ_max)+Fint(l1+l2,l3,δ,τ_max))

end

function Lint(l1, l2, δ, τ_max)
        
    L = 1/l1*(Gint(l2,δ,τ_max)-Fint(l1,l2,δ,τ_max))

end

function calcM0(x, prms::params{R,I}) where {R<:Real, I<:Integer}

    @unpack  δ_th, α, γ, W_bar, θ, θ0, κ_r, κ_β, σ2_r, σ2_β, τ_max = prms
    
    Z      = σ2_r*(γ/W_bar)*θ*Hint(κ_r,δ_th[1],δ_th[2],τ_max)/(1 + σ2_r*(γ/W_bar)*α*(Fint(κ_r,κ_r,δ_th[1],τ_max)-Fint(κ_r,x[1],δ_th[1],τ_max))/(x[1]-κ_r))
    M0obj = x[1] - κ_β + σ2_β*(γ/W_bar)*θ*Z*(Hint(κ_r,δ_th[1],δ_th[2],τ_max) - Hint(x[1],δ_th[1],δ_th[2],τ_max) - Z*(Fint(κ_r,κ_r,δ_th[1],τ_max) + Fint(x[1],x[1],δ_th[1],τ_max) - 2*Fint(κ_r,x[1],δ_th[1],τ_max))/(x[1] - κ_r))/(x[1] - κ_r)
    
    return M0obj

end
