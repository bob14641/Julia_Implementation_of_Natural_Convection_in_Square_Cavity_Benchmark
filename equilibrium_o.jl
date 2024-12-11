function equilibrium_o(nnx,nny,nnode,ui,uj,temp,rho,ndir,vp,wi,invcs2,invcs4_2,invcs2_2,f1,t1,f2,t2)

    for x = 1 : nnx
        for y = 1 : nny

            in = indicies_scalar(nnx,x,y)

            uilc = ui[in]
            ujlc = uj[in]
            rholc = rho[in]
            tlc = temp[in]
            udotu = (uilc * uilc) + (ujlc * ujlc)

            for q = 1 : ndir
                ina = indicies_field(nnx,nnode,x,y,q)
                dirv = (vp[q,1] * uilc) + (vp[q,2] * ujlc)
                feq0 = (wi[q] * rholc) * (1.0 + (invcs2 * dirv) + (invcs4_2 * (dirv * dirv)) - invcs2_2 * udotu)
                f1[ina] = feq0
                f2[ina] = feq0
                teq0 = ((wi[q] * tlc) * (1.0 + (invcs2 * dirv) + (invcs4_2 * (dirv * dirv)) - invcs2_2 * udotu))
                t1[ina] = teq0
                t2[ina] = teq0
            end
        end
    end
end