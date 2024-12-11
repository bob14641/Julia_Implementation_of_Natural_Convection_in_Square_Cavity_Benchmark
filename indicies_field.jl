function indicies_field(nnx,nnode,x,y,ndir)

    fidx = ((nnx * 9) * (x - 1) + ((y - 1) * 9) + ndir)

    return fidx
    
end