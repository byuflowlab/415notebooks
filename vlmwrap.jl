import VLM
using PyPlot

function simpleversion(AR, λ, Λ, θr, θt, α)

    b = 1.0  # arbitrary constant
    ϕ = 0.0
    npanels = 100
    duplicate = false
    spacing = "uniform"
    panels = VLM.simplewing(b, AR, λ, Λ, ϕ, θr, θt, npanels, duplicate, spacing)

    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    vother = nothing
    fs = VLM.Freestream(α, beta, Omega, vother)

    Sref = b^2/AR
    cref = Sref/b
    bref = b
    rcg = [0.0, 0.0, 0.0]
    ref = VLM.Reference(Sref, cref, bref, rcg)

    symmetric = true
    CF, CM, ymid, zmid, l, cl, dCF, dCM = VLM.run(panels, ref, fs, symmetric)
    CDi, CY, CL = CF
    Cl, Cm, Cn = CM

    einv = CL^2/(pi*AR*CDi)
    rootx = panels[1].rl[1] + 3.0/4*panels[1].chord
    tipx = panels[end].rr[1] + 3.0/4*panels[end].chord
    maxx = max(rootx, tipx)

    subplot(121)
    VLM.visualizeoutline2d(panels)
    text(-b/2.0, -maxx - 0.1, @sprintf("CL = %0.3f", CL))
    text(-b/2.0, -maxx - 0.2, @sprintf("CDi_inv = %0.5f", CDi))
    text(-b/2.0, -maxx - 0.3, @sprintf("e_inv = %0.3f", einv))

    subplot(122)
    plot(ymid/b, l)
    plot(ymid/b, cl)
    # plot([0, 0.5], [clmax, clmax], "k--")
    ylim([0, 1.3])
    xlabel("y/b")
    legend(["lift", "cl"])

end

# AR = 8.0
# λ = 0.55
# Λ = 25*pi/180
# θr = 0.0*pi/180
# θt = 0.0*pi/180
# α = 10*pi/180
# figure(figsize=(8, 3))
# simpleversion(AR, λ, Λ, θr, θt, α)   
# gcf()


function definewing(xle, yle, zle, chord, theta, npanels)

    bref = 2*(yle[end] - yle[1])
    
    n = length(xle)
    Sref = 0.0
    cref = 0.0  # use cmac
    for i = 1:n-1
        ds = sqrt((yle[i+1]-yle[i])^2 + (zle[i+1]-zle[i])^2)
        cr = chord[i]
        ct = chord[i+1]
        cmaci = 2.0/3 * (cr + ct - cr*ct/(cr + ct))
        Si = 0.5*(cr + ct) * ds
        cref += Si*cmaci
        Sref += Si
    end
    cref /= Sref
    Sref *= 2

    rcg = [xle[1] + 0.25*chord[1], yle[1], zle[1]]
    
    duplicate = false
    spacing = "uniform"
    wing = VLM.linearsections(xle, yle, zle, chord, theta, npanels, duplicate, spacing)

    return wing, Sref, cref, bref, rcg
end

function analyzewing(wing, Sref, cref, bref, rcg, alpha, clmax)

    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    vother = nothing
    fs = VLM.Freestream(alpha, beta, Omega, vother)

    
    ref = VLM.Reference(Sref, cref, bref, rcg)

    symmetric = true
    CF, CM, ymid, zmid, l, cl, dCF, dCM = VLM.run(wing, ref, fs, symmetric)
    CDi, CY, CL = CF
    Cl, Cm, Cn = CM
    AR = bref^2/Sref
    einv = CL^2/(pi*AR*CDi)

    figure(figsize=(8, 3))
    subplot(121)
    VLM.twoview(wing)
    
    subplot(122)
    plot(ymid/bref, cl)
    plot(ymid/bref, l)
    plot([0, 0.5], [clmax, clmax], "k--")
    xlabel("y/b")
    legend(["lift", "cl"])

    println("")
    println("")
    # println("CL = ", CL)
    @printf("CL = %0.3f\n", CL)
    @printf("CDi_inv = %0.5f\n", CDi)
    @printf("e_inv = %0.3f\n", einv)
    @printf("Cm = %0.3f\n", Cm)
    @printf("S_ref = %0.3f\n", Sref)
    @printf("c_ref = %0.3f\n", cref)


end


function wing(xle, yle, zle, chord, theta, npanels, alpha, clmax)

    wing, Sref, cref, bref, rcg = definewing(xle, yle, zle, chord, theta, npanels)
    analyzewing(wing, Sref, cref, bref, rcg, alpha, clmax)
    
end



function biplane(xle, yle, zle, chord, theta, npanels, voffset, alpha, clmax)

    wing1, Sref, cref, bref, rcg = definewing(xle, yle, zle, chord, theta, npanels)
    wing2, Sref, cref, bref, rcg = definewing(xle, yle, zle + voffset, chord, theta, npanels)
    wing = [wing1; wing2]
    analyzewing(wing, Sref, cref, bref, rcg, alpha, clmax)
end


# xle = [0.0; 0.4]
# yle = [0.0; 7.5]
# zle = [0.0; 0.0]
# chord = [2.2; 1.8]
# theta = [2.0*pi/180; 2.0*pi/180]
# npanels = [11]
# alpha = 9.0*pi/180
# clmax = 1.3

# wing(xle, yle, zle, chord, theta, npanels, alpha, clmax)
# gcf()
