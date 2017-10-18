using PyPlot
using Interpolations
using Roots

struct motordef
    KvRPM::Float64
    i0::Float64
    R::Float64
    vbatt::Float64
end

struct propdef
    data::Array{Float64, 2}
    D::Float64
    rho::Float64
end

struct acdef
    CDp::Float64
    einv::Float64
    Sref::Float64
    b::Float64
    nmotors
    mass::Float64
end

function interp1(xpt, ypt, x)
    """helper function"""
    intf = interpolate((xpt,), ypt, Gridded(Linear()))
    return intf[x]
end


function motor(m::motordef, OmegaRPM, throttle)
    """units: RPM, Amp, Ohm, Voltage, RPM"""

    # unit conversion
    Kv = m.KvRPM*pi/30
    Omega = OmegaRPM*pi/30

    v = m.vbatt*throttle  # assumes linear
    i = (v - Omega/Kv)/m.R  # current
    Q = (i - m.i0)/Kv  # torque
    Pout = Q*Omega  # power out
    eta = Pout/(i*v)  # efficiency

    return i, Q, Pout, eta
end


function prop(p::propdef, OmegaRPM, V)

    # unit conversion
    Omega = OmegaRPM*pi/30
    n = Omega/(2*pi)

    # extract data
    J_data = p.data[:, 1]
    CT_data = p.data[:, 2]
    CP_data = p.data[:, 3]
    eta_data = p.data[:, 4]

    # interpolate off of charts
    J = V/(n*p.D)  # advance ratio
    CT = interp1(J_data, CT_data, J)
    CP = interp1(J_data, CP_data, J)
    eta = interp1(J_data, eta_data, J)

    CQ = CP/(2*pi)

    T = CT * p.rho .* n.^2 * p.D^4
    Q = CQ * p.rho .* n.^2 * p.D^5
    # Pout = CP * rho .* n.^3 * D^5
    # Pout = T*V

    return Q, T, eta
end


function motorpropsweep(m::motordef, p::propdef, Omegamin, Omegamax, throttle, V)

    # evalute at a range of rotation speeds
    no = 100
    Omegavec = linspace(Omegamin, Omegamax, no)  # RPM
    Qmvec = zeros(no)
    etamvec = zeros(no)
    Qpvec = zeros(no)
    Tpvec = zeros(no)
    etapvec = zeros(no)

    # run a sweep
    for i = 1:no
        im, Qmvec[i], Poutm, etamvec[i] = motor(m, Omegavec[i], throttle)
        Qpvec[i], Tpvec[i], etapvec[i] = prop(p, Omegavec[i], V)
    end

    # solve for when torque's are equal
    function difference(OmegaRPM)
        im, Qm, Poutm, etam = motor(m, OmegaRPM, throttle)
        Qp, Tp, etap = prop(p, OmegaRPM, V)
        return Qm-Qp
    end

    OmegaRPMO = fzero(difference, 1.0, 100000)

    # evaluate models at this Omega
    imO, QmO, PoutmO, etamO = motor(m, OmegaRPMO, throttle)
    QpO, TpO, etapO = prop(p, OmegaRPMO, V)

    # total efficiency
    eta = etamO*etapO

    figure()
    plot(Omegavec, Qmvec*1e3)
    plot(Omegavec, Qpvec*1e3)
    plot([OmegaRPMO; OmegaRPMO], [0; maximum(Qmvec*1e3)], "k--")
    xlabel("Omega (RPM)")
    ylabel("Q (N-mm)")
    ylim([0, maximum(Qmvec*1e3)])
    legend(["motor", "prop"])

    figure()
    plot(Omegavec, etamvec)
    plot(Omegavec, etapvec)
    plot([OmegaRPMO; OmegaRPMO], [0; 1.0], "k--")
    xlabel("Omega (RPM)")
    ylabel(L"\eta")
    ylim([0, 1])
    legend(["motor", "prop"])

    figure()
    plot(Omegavec, Tpvec)
    plot([OmegaRPMO; OmegaRPMO], [0; maximum(Tpvec)], "k--")
    xlabel("Omega (RPM)")
    ylabel("T (N)")

    return eta, TpO

end


function steadylevel(m::motordef, p::propdef, ac::acdef, V)

    # compute drag (and thus required thrust)
    mass = ac.mass/1e3  # convert to kg
    g = 9.81
    W = mass*g
    q = 0.5*p.rho*V^2  # dynamic pressure
    AR = ac.b^2/ac.Sref  # aspect ratio
    e = 1.0 / (1.0/ac.einv + 0.38*ac.CDp*pi*AR)  # Oswald efficiency factor
    D = ac.CDp*q*ac.Sref + W^2/(q*pi*ac.b^2*e)  # total drag
    Tdesired = D/ac.nmotors  # desired thrust per motor

    # find location where thrust = Tdesired
    function difference(OmegaRPM)
        Q, T, eta = prop(p, OmegaRPM, V)
        return T - Tdesired
    end

    # find location where thrust is thrust-desired
    OmegaRPM = fzero(difference, 1.0, 100000)

    # evaluate propeller torque
    Qp, Tp, etap = prop(p, OmegaRPM, V)

    # unit conversions
    Kv = m.KvRPM*pi/30
    Omega = OmegaRPM*pi/30

    # requisite voltage and throttle
    v = (Kv*Qp + m.i0)*m.R + Omega/Kv
    throttle = v/m.vbatt

    # evaluate motor
    im, Qm, Poutm, etam = motor(m, OmegaRPM, throttle)

    # Pelec = Pm/etam
    # Pthrust = T*V

    return etap, etam, throttle, OmegaRPM, im, Tp
end


function steadylevelsweep(m::motordef, p::propdef, ac::acdef, Vmin, Vmax, Omegamax, imax, Tmax)

    nv = 50
    Vvec = linspace(Vmin, Vmax, nv)
    throttle = zeros(nv)
    etap = zeros(nv)
    etam = zeros(nv)
    OmegaRPM = zeros(nv)
    current = zeros(nv)
    Tp = zeros(nv)

    for i = 1:nv
        etap[i], etam[i], throttle[i], OmegaRPM[i], current[i], Tp[i] = steadylevel(m, p, ac, Vvec[i])
    end

    figure()
    plot(Vvec, throttle)
    xlabel("flight speed")
    ylabel("throttle")
    ylim([0, 1])

    figure()
    plot(Vvec, etam)
    plot(Vvec, etap)
    plot(Vvec, etam.*etap)
    xlabel("flight speed")
    ylabel("efficiency")
    legend(["motor", "prop", "combined"])
    ylim([0, 1])


    figure()
    plot(Vvec, OmegaRPM)
    xlabel("flight speed")
    ylabel("Omega")
    ylim([0, Omegamax])

    figure()
    plot(Vvec, current)
    xlabel("flight speed")
    ylabel("motor current (amps)")
    ylim([0, imax])

    figure()
    plot(Vvec, Tp)
    xlabel("flight speed")
    ylabel("Thrust per motor (N)")
    ylim([0, Tmax])

end
