"""
    window_joglekar(w, w_on, w_off, p=1)

Joglekar window function to bound state variable.
f(w) = 1 - ((2w - w_on - w_off)/(w_off - w_on))^(2p)
"""
function window_joglekar(w::Real, w_on::Real, w_off::Real, p::Int=1)
    w_norm = (2w - w_on - w_off) / (w_off - w_on)
    return 1 - w_norm^(2p)
end

"""
    window_biolek(w, w_on, w_off, dw, p=1)

Biolek window function (direction-dependent).
"""
function window_biolek(w::Real, w_on::Real, w_off::Real, dw::Real, p::Int=1)
    if dw >= 0
        return 1 - ((w - w_on) / (w_off - w_on))^(2p)
    else
        return 1 - ((w_off - w) / (w_off - w_on))^(2p)
    end
end
