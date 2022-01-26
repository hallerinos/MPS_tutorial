function tensor_to_matrix(tensor::ITensor)::Matrix
    indices = inds(tensor)
    primed_indices = indices[1:2:end]
    normal_indices = indices[2:2:end]
    arr = Array(tensor, [primed_indices..., normal_indices...]...)
    mat = reshape(arr, prod(dims(primed_indices)), prod(dims(normal_indices)))
    return mat
end