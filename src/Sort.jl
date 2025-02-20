
@propagate_inbounds function bubbleSort!(a::AbstractVector; lt = isless)
    n = length(a)
    for i = 1:(n - 1)
        swapped = false
        for j = 1:(n - i)
            if lt(a[j + 1], a[j])
                a[j], a[j + 1] = (a[j + 1], a[j])
                swapped = true
            end
        end

        # If no two elements were swapped
        # by inner loop, then break
        if !swapped
            break
        end
    end
end

@propagate_inbounds function bubbleSortPerm!(
    ix::AbstractVector,
    a::AbstractVector;
    lt = isless,
    initialized = false,
)
    if !initialized
        ix .= LinearIndices(a)
    end
    n = length(a)
    for i = 1:(n - 1)
        swapped = false
        for j = 1:(n - i)
            if lt(a[ix[j + 1]], a[ix[j]])
                ix[j], ix[j + 1] = (ix[j + 1], ix[j])
                swapped = true
            end
        end

        # If no two elements were swapped
        # by inner loop, then break
        if !swapped
            break
        end
    end
end

@propagate_inbounds function cocktailSort!(a::AbstractVector; lt = isless)
    n = length(a)
    swapped = true
    start = 1
    stop = n - 1

    while swapped
        # reset the swapped flag on entering
        # the loop, because it might be true from
        # a previous iteration.
        swapped = false

        # loop from left to right same as the bubble sort
        for i = start:stop
            if lt(a[i + 1], a[i])
                a[i], a[i + 1] = (a[i + 1], a[i])
                swapped = true
            end
        end

        # if nothing moved, then array is sorted.
        if !swapped
            break
        end

        # otherwise, reset the swapped flag so that it
        # can be used in the next stage
        swapped = false

        # move the end point back by one, because
        # item at the end is in its rightful spot
        stop -= 1

        # from right to left, doing the
        # same comparison as in the previous stage
        for i = stop:-1:start
            if lt(a[i + 1], a[i])
                a[i], a[i + 1] = (a[i + 1], a[i])
                swapped = true
            end
        end

        # increase the starting point, because
        # the last stage would have moved the next
        # smallest number to its rightful spot.
        start += 1
    end
end


# To find gap between elements
@inline function getNextGap(gap::Int)
    # Shrink gap by Shrink factor
    gap = unsafe_trunc(Int, (gap * 10) / 13)

    if gap < 1
        return 1
    end
    return gap
end

@propagate_inbounds function combSort!(a::AbstractVector; lt = isless)
    # Initialize gap
    n = length(a)
    gap = n

    # Initialize swapped as true to make sure that loop runs
    swapped = true

    # Keep running while gap is more than 1 and last
    # iteration caused a swap
    while gap != 1 || swapped == true
        # Find next gap
        gap = getNextGap(gap)

        # Initialize swapped as false so that we can
        # check if swap happened or not
        swapped = false

        # Compare all elements with current gap
        for i = 1:(n - gap)
            if lt(a[i + gap], a[i])
                a[i], a[i + gap] = (a[i + gap], a[i])
                swapped = true
            end
        end
    end
end
