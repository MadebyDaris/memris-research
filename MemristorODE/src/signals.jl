"""
    triangular_wave(t; amplitude, period)

Generate triangular wave voltage signal.
"""
function triangular_wave(t; amplitude=1.0, period=1e-3)
    phase = mod(t, period) / period
    if phase < 0.25
        return 4 * amplitude * phase
    elseif phase < 0.75
        return amplitude * (2 - 4 * phase)
    else
        return 4 * amplitude * (phase - 1)
    end
end

"""
    sinusoidal_wave(t; amplitude, frequency)

Generate sinusoidal voltage signal.
"""
function sinusoidal_wave(t; amplitude=1.0, frequency=1e3)
    return amplitude * sin(2π * frequency * t)
end

"""
    pulse_train(t; V_high, V_low, t_high, t_low, t_rise=0.0)

Generate pulse train for SET/RESET programming.
"""
function pulse_train(t; V_high=1.0, V_low=-1.0, t_high=1e-6, t_low=1e-6, t_rise=0.0)
    T = t_high + t_low + 2*t_rise
    phase = mod(t, T)

    if phase < t_rise
        return V_low + (V_high - V_low) * (phase / t_rise)
    elseif phase < t_rise + t_high
        return V_high
    elseif phase < 2*t_rise + t_high
        return V_high + (V_low - V_high) * ((phase - t_rise - t_high) / t_rise)
    else
        return V_low
    end
end
