import Adapt

function binEvents!(h::Hist3, events::AbstractArray, transforms::Array1{SquareMatrix3c})
    JACC.parallel_for(
        size(events, 2),
        (i, t) -> begin
            @inbounds begin
                for n in 1:t.len
                    op = t.transforms[n]
                    v = op * C3[t.events[6, i], t.events[7, i], t.events[8, i]]
                    atomic_push!(t.h, v[1], v[2], v[3], t.events[1, i])
		end
            end
        end,
        (h, events, transforms, len = length(transforms)),
    )
end
