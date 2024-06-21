
@propagate_inbounds function bubbleSort!(a::AbstractVector; lt = isless, iPerm = nothing)
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
                # a[j], a[j + 1] = (a[j + 1], a[j])
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

# void bubbleSort(int arr[], int n)
# {
#     int i, j;
#     bool swapped;
#     for (i = 0; i < n - 1; i++) {
#         swapped = false;
#         for (j = 0; j < n - i - 1; j++) {
#             if (arr[j] > arr[j + 1]) {
#                 swap(arr[j], arr[j + 1]);
#                 swapped = true;
#             }
#         }

#         // If no two elements were swapped
#         // by inner loop, then break
#         if (swapped == false)
#             break;
#     }
# }
