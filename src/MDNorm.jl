import Base: @propagate_inbounds

mutable struct MDNorm
    hX::Array1r
    kX::Array1r
    lX::Array1r
    eX::Array1r
    origin::Crd3
    boxLength::Crd3
    indexMax::Id3
    indexMaker::Id3

    function MDNorm(hx::Array1r, kx::Array1r, lx::Array1r)
        indexMax = I3[200, 200, 1]
        new(
            hx,
            kx,
            lx,
            Array1r(),
            C3[-10.0, -10.0, -0.1],
            C3[0.1, 0.1, 0.2],
            indexMax,
            setUpIndexMaker(indexMax),
        )
    end
end

function MDNorm()
    MDNorm(Array1r(), Array1r(), Array1r())
end

function MDNorm(hx::Vector, kx::Vector, lx::Vector)
    MDNorm(Array1r(hx), Array1r(kx), Array1r(lx))
end

"""
    calculateIntersections!() -> Vector{Vec4}

Calculate the points of intersection for the given detector with cuboid
surrounding the detector position in HKL.

# Arguments
- `histogram`: umm...
- `theta::Float64`: Polar angle with detector
- `phi::Float64`: Azimuthal angle with detector
- `transform::SquareMatrix3`: Matrix to convert frm Q_lab to HKL ``(2Pi*R *UB*W*SO)^{-1}``
- `lowvalue::Float64`: The lowest momentum or energy transfer for the trajectory
- `highvalue::Float64`: The highest momentum or energy transfer for the trajectory
"""
@propagate_inbounds function calculateIntersections!(
    mdn::MDNorm,
    histogram::THistogram,
    theta::Float64,
    phi::Float64,
    transform::SquareMatrix3,
    lowvalue::Float64,
    highvalue::Float64,
    intersections::Vector{Vec4},
) where {THistogram}
    sin_theta = sin(theta)
    qout = V3[sin_theta * cos(phi), sin_theta * sin(phi), cos(theta)]
    qin = V3[0, 0, 1]

    qout = transform * qout
    qin = transform * qin

    kimin = lowvalue
    kimax = highvalue
    kfmin = kimin
    kfmax = kimax

    # hx = Vector(mdn.hX)
    # kx = Vector(mdn.kX)
    # lx = Vector(mdn.lX)
    hx = mdn.hX
    kx = mdn.kX
    lx = mdn.lX
    hNPts = length(hx)
    kNPts = length(kx)
    lNPts = length(lx)

    hStart = qin[1] * kimin - qout[1] * kfmin
    hEnd = qin[1] * kimax - qout[1] * kfmax
    hStartIdx = binindex(histogram, hStart, 0, 0)[1]
    hEndIdx = binindex(histogram, hEnd, 0, 0)[1]
    hStartIdx = clamp(hStartIdx, 0, hNPts)
    hEndIdx = clamp(hEndIdx, 0, hNPts)
    if hStartIdx > hEndIdx
        hStartIdx, hEndIdx = swap(hStartIdx, hEndIdx)
    end
    hStartIdx += 1

    kStart = qin[2] * kimin - qout[2] * kfmin
    kEnd = qin[2] * kimax - qout[2] * kfmax
    kStartIdx = binindex(histogram, 0, kStart, 0)[2]
    kEndIdx = binindex(histogram, 0, kEnd, 0)[2]
    kStartIdx = clamp(kStartIdx, 0, kNPts)
    kEndIdx = clamp(kEndIdx, 0, kNPts)
    if kStartIdx > kEndIdx
        kStartIdx, kEndIdx = swap(kStartIdx, kEndIdx)
    end
    kStartIdx += 1

    lStart = qin[3] * kimin - qout[3] * kfmin
    lEnd = qin[3] * kimax - qout[3] * kfmax
    lStartIdx = binindex(histogram, 0, 0, lStart)[3]
    lEndIdx = binindex(histogram, 0, 0, lEnd)[3]
    lStartIdx = clamp(lStartIdx, 0, lNPts)
    lEndIdx = clamp(lEndIdx, 0, lNPts)
    if lStartIdx > lEndIdx
        lStartIdx, lEndIdx = swap(lStartIdx, lEndIdx)
    end
    lStartIdx += 1

    empty!(intersections)
    sizehint!(intersections, hNPts + kNPts + lNPts + 2)

    # calculate intersections with planes perpendicular to h
    fmom = (kfmax - kfmin) / (hEnd - hStart)
    fk = (kEnd - kStart) / (hEnd - hStart)
    fl = (lEnd - lStart) / (hEnd - hStart)
    for i = hStartIdx:hEndIdx
        hi = hx[i]
        # if hi is between hStart and hEnd, then ki and li will be between
        # kStart, kEnd and lStart, lEnd and momi will be between kfmin and
        # kfmax
        ki = fk * (hi - hStart) + kStart
        li = fl * (hi - hStart) + lStart
        if (ki >= kx[1]) && (ki <= kx[kNPts]) && (li >= lx[1]) && (li <= lx[lNPts])
            momi = fmom * (hi - hStart) + kfmin
            push!(intersections, V4[hi, ki, li, momi])
        end
    end

    # calculate intersections with planes perpendicular to k
    fmom = (kfmax - kfmin) / (kEnd - kStart)
    fh = (hEnd - hStart) / (kEnd - kStart)
    fl = (lEnd - lStart) / (kEnd - kStart)
    for i = kStartIdx:kEndIdx
        ki = kx[i]
        # if ki is between kStart and kEnd, then hi and li will be between
        # hStart, hEnd and lStart, lEnd and momi will be between kfmin and
        # kfmax
        hi = fh * (ki - kStart) + hStart
        li = fl * (ki - kStart) + lStart
        if (hi >= hx[1]) && (hi <= hx[hNPts]) && (li >= lx[1]) && (li <= lx[lNPts])
            momi = fmom * (ki - kStart) + kfmin
            push!(intersections, V4[hi, ki, li, momi])
        end
    end

    # calculate intersections with planes perpendicular to l
    fmom = (kfmax - kfmin) / (lEnd - lStart)
    fh = (hEnd - hStart) / (lEnd - lStart)
    fk = (kEnd - kStart) / (lEnd - lStart)

    for i = lStartIdx:lEndIdx
        li = lx[i]
        hi = fh * (li - lStart) + hStart
        ki = fk * (li - lStart) + kStart
        if (hi >= hx[1]) && (hi <= hx[hNPts]) && (ki >= kx[1]) && (ki <= kx[kNPts])
            momi = fmom * (li - lStart) + kfmin
            push!(intersections, V4[hi, ki, li, momi])
        end
    end

    # endpoints
    if (hStart >= hx[1]) &&
       (hStart <= hx[hNPts]) &&
       (kStart >= kx[1]) &&
       (kStart <= kx[kNPts]) &&
       (lStart >= lx[1]) &&
       (lStart <= lx[lNPts])
        push!(intersections, V4[hStart, kStart, lStart, kfmin])
    end
    if (hEnd >= hx[1]) &&
       (hEnd <= hx[hNPts]) &&
       (kEnd >= kx[1]) &&
       (kEnd <= kx[kNPts]) &&
       (lEnd >= lx[1]) &&
       (lEnd <= lx[lNPts])
        push!(intersections, V4[hEnd, kEnd, lEnd, kfmax])
    end

    # sort intersections by final momentum
    sort!(intersections, lt = (v1, v2) -> v1[4] < v2[4])

    return intersections
end


"""
Calculate the normalization among intersections on a single detector in 1
specific SpectrumInfo/ExperimentInfo.

# Arguments
- `intersections::Vector{Vec4}`: intersections
- `solid::ScalarType`: proton charge
- `yValues::Vector{ScalarType}`: diffraction intersection integral and common to sample and background
"""
@propagate_inbounds function calculateSingleDetectorNorm!(
    mdn::MDNorm,
    intersections::Vector{Vec4},
    solid::ScalarType,
    yValues::Vector{ScalarType},
    histogram::THistogram,
) where {THistogram}
    for i = 2:length(intersections)
        curX = intersections[i]
        prevX = intersections[i - 1]

        # The full vector isn't used so compute only what is necessary
        # If the difference between 2 adjacent intersection is trivial, no
        # intersection normalization is to be calculated
        # diffraction
        delta = curX[4] - prevX[4]
        eps = 1.0e-7
        if delta < eps
            continue # Assume zero contribution if difference is small
        end

        # Average between two intersections for final position
        # [Task 89] Sample and background have same 'pos[]'
        pos = V3[
            0.5 * (curX[1] + prevX[1]),
            0.5 * (curX[2] + prevX[2]),
            0.5 * (curX[3] + prevX[3]),
        ]

        # Diffraction
        # signal = integral between two consecutive intersections
        signal = (yValues[i] - yValues[i - 1]) * solid

        atomic_push!(histogram, pos[1], pos[2], pos[3], signal)
    end
end

struct XValues
    ix::Vector{Vec4}
end

@propagate_inbounds @inline Base.getindex(xValues::XValues, i::Integer) = xValues.ix[i][4]
@inline Base.length(xValues::XValues) = length(xValues.ix)

"""
Linearly interpolate between the points in integrFlux at xValues and save the
results in yValues.

# Arguments
- `xValues`: X-values at which to interpolate
- `integrFlux`: A workspace with the spectra to interpolate
- `sp::SizeType`: A workspace index for a spectrum in integrFlux to interpolate.
"""
@propagate_inbounds function calculateIntegralsForIntersections(
    xValues::XValues,
    integrFlux_x::AbstractRange,
    integrFlux_y::Vector{Vector{ScalarType}},
    sp::SizeType,
)

    # the x-data from the workspace
    xData = integrFlux_x
    xStart = first(xData)
    xEnd = last(xData)

    # the values in integrFlux are expected to be integrals of a non-negative
    # function
    # ie they must make a non-decreasing function
    yData = integrFlux_y[sp]
    yMin = 0.0
    yMax = last(yData)

    nData = length(xValues)

    # all integrals below xStart must be 0
    if xValues[nData] < xStart
        return fill(yMin, nData)
    end

    # all integrals above xEnd must be equal to yMax
    if xValues[1] > xEnd
        return fill(yMax, nData)
    end

    yValues = Vector{ScalarType}(undef, nData)

    i = 1
    # integrals below xStart must be 0
    while i < nData && xValues[i] <= xStart
        yValues[i] = yMin
        i += 1
    end

    iMax = nData
    # integrals above xEnd must be yMax
    while iMax > i && xValues[iMax] >= xEnd
        yValues[iMax] = yMax
        iMax -= 1
    end

    j = 2
    while i <= iMax
        xi = xValues[i]
        while xi > xData[j]
            j += 1
        end
        # interpolate between the consecutive points
        x0 = xData[j - 1]
        x1 = xData[j]
        y0 = yData[j - 1]
        y1 = yData[j]
        yValues[i] = lerp(y0, y1, (xi - x0) / (x1 - x0))
        i += 1
    end

    return yValues
end

"""
Calculate the diffraction MDE's intersection integral of a certain
detector/spectru

# Arguments
- `intersections::Vector{Vec4}`: vector of intersections
- `integrFlux`: integral flux workspace
- `wsIdx::SizeType`: workspace index
"""
@propagate_inbounds function calculateDiffractionIntersectionIntegral(
    intersections::Vector{Vec4},
    integrFlux_x::AbstractRange,
    integrFlux_y::Vector{Vector{ScalarType}},
    wsIdx::SizeType,
)
    # xValues = map(inte -> inte[4], intersections)
    xValues = XValues(intersections)
    yValues = calculateIntegralsForIntersections(xValues, integrFlux_x, integrFlux_y, wsIdx)
    return (xValues, yValues)
end
