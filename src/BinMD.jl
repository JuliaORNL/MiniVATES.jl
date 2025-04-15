import Adapt

function binEvents!(h::Hist3, events::AbstractArray, transforms::Array2c)
    JACC.parallel_for(
    size(events, 1),
        (i, t) -> begin
            transforms_shared = JACC.shared(t.transforms)
            @inbounds begin
		for n in 1:t.hjkl
	            tmp = (n-1)*3
                    v = StaticArrays.MVector{3,CoordType}(0, 0, 0)
		    for j in 1:3
			tmp2 = j * t.asdf;
                        for k in 1:3
                            v[j] += t.events[i, k] * transforms_shared[tmp + k + tmp2]
                        end
                    end
                    atomic_push!(t.h, v[1], v[2], v[3], 1.)
                end
            end
        end,
        (h = h, events, transforms, asdf = size(transforms, 1), hjkl= div(size(transforms,1),3)),
    )
end

