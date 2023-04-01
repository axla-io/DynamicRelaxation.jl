using Base.Iterators

# Define the ranges
range1 = 1:4:10
range2 = 2:4:10

# Define a function that interleaves the elements from two iterators
function interleave(itr1, itr2)
    # Define a nested generator function that generates the interleaved elements
    function interleaved_gen()
        # Create two iterators from the input iterators
        itr1_copy = itr1
        itr2_copy = itr2
        
        # Loop until both iterators are exhausted
        while (itr1_copy !== nothing) || (itr2_copy !== nothing)
            # Yield the next element from the first iterator, if it exists
            if itr1_copy !== nothing
                Base.Iterators.yield(first(itr1_copy))
                itr1_copy = next(itr1_copy, nothing)
            end
            
            # Yield the next element from the second iterator, if it exists
            if itr2_copy !== nothing
                Base.Iterators.yield(first(itr2_copy))
                itr2_copy = next(itr2_copy, nothing)
            end
        end
    end
    
    # Return an iterator that uses the nested generator function
    return Stateful(interleaved_gen())
end

# Interleave the elements from the two ranges
interleaved_range = interleave(range1, range2)

# Test the iterator
for i in interleaved_range
    println(i)
end