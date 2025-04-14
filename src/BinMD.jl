import Adapt

function binEvents!(h::Hist3, events::AbstractArray, transforms::Array1{SquareMatrix3c})
    JACC.parallel_for(
        (length(transforms), size(events, 2)),
        (n, i, t) -> begin
            @inbounds begin
                op = t.transforms[n]
                v = op * C3[t.events[i, 1], t.events[i, 2], t.events[i, 3]]
                atomic_push!(t.h, v[1], v[2], v[3], t.events[1, i])
            end
        end,
        (h = h, events, transforms),
    )
end
