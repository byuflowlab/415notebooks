function vortex(xA, yA, xB, yB, xC, yC)

    m = length(yC)
    n = length(yA)

    w = zeros(m,n)

    for i = 1:m
        for j = 1:n
            # first define some repeatedly used values
            x1 = xC[i] - xA[j]
            x2 = xC[i] - xB[j]
            x21 = xB[j] - xA[j]
            y1 = yC[i] - yA[j]
            y2 = yC[i] - yB[j]
            y21 = yB[j] - yA[j]

            denom1 = (x1*y2-x2*y1)^2
            k1 = x1*y2 - x2*y1
            frac1_1 = (x21*x1 + y21*y1)/sqrt(x1^2 + y1^2)
            frac1_2 = (x21*x2 + y21*y2)/sqrt(x2^2 + y2^2)
            frac1 = frac1_1 - frac1_2
            w[i,j] = 1.0/(4*pi)*k1/denom1*frac1
        end
    end

    return w
end

function getAIC(QCx, QCy, CPx, CPy)

    # first calculate induced velocities from vorticies on right side of wing.
    m = length(CPy)
    n = length(QCy)-1

    AIC = zeros(m, n)


    # Induced velocities due to bound vortex (BC)
    AIC += vortex(QCx[1:n], QCy[1:n], QCx[2:n+1], QCy[2:n+1], CPx, CPy)

    # Induced velocities due to trailing vortex from A
    xINF = QCx + 1e9
    yINF = QCy

    AIC += vortex(xINF[1:n], yINF[1:n], QCx[1:n], QCy[1:n], CPx, CPy)

    # Induced velocities due to trailing vortex from D
    AIC += vortex(QCx[2:n+1],QCy[2:n+1],xINF[2:n+1],yINF[2:n+1],CPx,CPy)


    # ---------- Now repeat with contribution from left side of the wing -----
    # y switched signs (except control points of course)
    QCy = -QCy
    yINF = -yINF
    # Also all induced velocities have opposite sign because vorticies travel
    # in opposite direction.

    # Induced velocities due to bound vortex (BC)
    AIC -= vortex(QCx[1:n],QCy[1:n],QCx[2:n+1],QCy[2:n+1],CPx,CPy)

    # Induced velocities due to trailing vortex from A

    AIC -= vortex(xINF[1:n],yINF[1:n],QCx[1:n],QCy[1:n],CPx,CPy)

    # Induced velocities due to trailing vortex from D
    AIC -= vortex(QCx[2:n+1],QCy[2:n+1],xINF[2:n+1],yINF[2:n+1],CPx,CPy)

    return AIC
end

function getDIC(yFF, rho)

    n = length(yFF)-1

    # ------------ define useful variables ------------------
    yc = 1/2*(yFF[2:n+1] + yFF[1:n]) # center of panels
    nz = yFF[2:n+1] - yFF[1:n]
    # ----------------------------------------------------

    DIC = zeros(n,n)
    #------- normal wash calculation ----------------
    for i = 1:n
        for j = 1:n
            ry = yFF[j] - yc[i]
            DIC[i, j] = rho/2/pi/ry*(-1)*nz[i]

            # add other side of panel
            ry = yFF[j+1] - yc[i]
            DIC[i, j] += rho/2/pi/ry*nz[i]
        end
    end

    yFF = -yFF # add left side

    for i = 1:n
        for j = 1:n
            ry = yFF[j] - yc[i]
            DIC[i, j] += rho/2/pi/ry*nz[i]

            # add other side of panel
            ry = yFF[j+1] - yc[i]
            DIC[i, j] += rho/2/pi/ry*(-1)*nz[i]
        end
    end

    return DIC
end

function geometry(chord, theta, b, chi, Lambda)

    N = 100
    N1 = round(chi*N)
    N2 = N - N1

    dy1 = chi*b/2
    dy2 = (1-chi)*b/2
    dx1 = dy1*tan(Lambda[1]*pi/180.0)
    dx2 = dy2*tan(Lambda[2]*pi/180.0)

    QCx1 = collect(0:N1)./N1 * dx1
    QCx2 = dx1 + collect(1:N2)./N2 * dx2
    QCx = [QCx1; QCx2]

    QCy1 = collect(0:N1)./N1 * dy1
    QCy2 = dy1 + collect(1:N2)./N2 * dy2
    QCy = [QCy1; QCy2]

    c1 = chord[1] + collect(0:N1)./N1 * (chord[2] - chord[1])
    c2 = chord[2] + collect(1:N2)./N2 * (chord[3] - chord[2])
    c = [c1; c2]
    ccp = 0.5*(c[1:N] + c[2:N+1])

    t1 = theta[1] + collect(0:N1)./N1 * (theta[2] - theta[1])
    t2 = theta[2] + collect(1:N2)./N2 * (theta[3] - theta[2])
    t = [t1; t2]
    thetacp = 0.5*(t[1:N] + t[2:N+1])

    CPy = 0.5*(QCy[1:N] + QCy[2:N+1])
    CPx = 0.5*(QCx[1:N] + QCx[2:N+1]) + 0.5*ccp

    ds = QCy[2:N+1] - QCy[1:N]

    return QCx, QCy, CPx, CPy, ccp, thetacp*pi/180.0, ds
end

function plotwing(chord, b, chi, Lambda)

    dy1 = chi*b/2
    dy2 = (1-chi)*b/2
    dx1 = dy1*tan(Lambda[1]*pi/180.0)
    dx2 = dy2*tan(Lambda[2]*pi/180.0)

    x1 = -chord[1]/4.0
    y1 = 0.0

    x2 = dx1 - chord[2]/4.0
    y2 = dy1

    x3 = dx1 + dx2 - chord[3]/4.0
    y3 = dy1 + dy2

    x4 = x3 + chord[3]
    y4 = y3

    x5 = x2 + chord[2]
    y5 = y2

    x6 = x1 + chord[1]
    y6 = y1

    plot([y1, y2, y3, y4, y5, y6], -[x1, x2, x3, x4, x5, x6], "k")
    plot(-[y1, y2, y3, y4, y5, y6], -[x1, x2, x3, x4, x5, x6], "k")
    axis("equal")
    axis("off")
end

function run(b, chi, chord, theta, Lambda, alpha, clmax)

    QCx, QCy, CPx, CPy, ccp, thetacp, ds = geometry(chord, theta, b, chi, Lambda)

    rho = 2.0
    V = 1.0
    q = 0.5*rho*V^2

    alpha *= pi/180
    Sref = 2*(0.5*(chord[1] + chord[2])*chi*b/2 + 0.5*(chord[2] + chord[3])*(1-chi)*b/2)
    AR = b^2/Sref

    AIC = getAIC(QCx, QCy, CPx, CPy)
    DIC = getDIC(QCy, rho)
    Vn = -V*(cos(alpha)*sin.(thetacp) + sin(alpha)*cos.(thetacp))
    gamma = AIC\Vn

    cl = rho*V*gamma./(q*ccp)
    l = cl.*ccp/mean(ccp)

    CL = 2*rho*V*sum(gamma.*ds)/(q*Sref)
    CDi = gamma'*DIC*gamma/(q*Sref)
    e = CL^2/(pi*AR*CDi)

    return CL, CDi, e, Sref, CPy/b, cl, l

end

function nondimrun(AR, tr, Lambda, twist, alpha, clmax)

    S = 1.0
    b = sqrt(S*AR)
    chi = 0.5
    cr = 2*S/(b*(1+tr))
    ct = cr*tr
    chord = [cr; 0.5*(cr + ct); ct]
    sweep = [Lambda; Lambda]

    CL, CDi, e, Sref, CPy, cl, l = run(b, chi, chord, twist, sweep, alpha, clmax)

    subplot(121)
    plotwing(chord, b, chi, sweep)
    text(-b/2.0, cr/4.0 + 0.6, @sprintf("CL = %0.3f", CL))
    text(-b/2.0, cr/4.0 + 0.4, @sprintf("CDi = %0.5f", CDi))
    text(-b/2.0, cr/4.0 + 0.2, @sprintf("e_inv = %0.3f", e))

    subplot(122)
    plot(CPy, l)
    plot(CPy, cl)
    xlabel("y/b")
    legend(["lift", "cl"])
end

function dimrun(b, chi, chord, theta, Lambda, alpha, clmax)
    CL, CDi, e, Sref, CPy, cl, l = run(b, chi, chord, theta, Lambda, alpha, clmax)

    figure(figsize=(8, 3))
    subplot(121)
    plotwing(chord, b, chi, Lambda)

    subplot(122)
    plot(CPy, l)
    plot(CPy, cl)
    plot([0, 0.5], [clmax, clmax], "k--")
    xlabel("y/b")
    legend(["lift", "cl", "clmax"])

    println("")
    println("")
    println("CL = ", CL)
    println("CDi = ", CDi)
    println("e_inv = ", e)
    println("S_ref = ", Sref)

end