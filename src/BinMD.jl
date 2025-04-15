import Adapt

function binEvents!(h::Hist3, events::AbstractArray, transforms::AbstractArray)
    JACC.parallel_for(
        size(events, 1),
        (i, t) -> begin
            transforms_shared = JACC.shared(t.transforms)
            @inbounds begin
		for n in 1:t.row_size
                    idx = (n-1)*9
                    v = StaticArrays.MVector{3,CoordType}(0, 0, 0)
		    for j in 1:3
                        for k in 1:3
			    v[k] += t.events[i, j] * transforms_shared[idx + k]
                        end
			idx += 3
                    end
                    atomic_push!(t.h, v[1], v[2], v[3], 1.)
                end
            end
        end,
        (h , events, transforms, row_size = size(transforms, 2)),
    )
end

