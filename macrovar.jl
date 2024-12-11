function macrovar(nnx,nny,nnode,f2,t2,rho,pres,cs2,fri,frj,ui,uj,temp,tmid,thot,tcold,betah,Fo,solidnode,uwxr,uwxl,uwxb,uwxt,uwyr,uwyl,uwyb,uwyt,dt)

    for x = 1 : nnx
        for y = 1 : nny
            in = indicies_scalar(nnx,x,y) 

            ina0 = indicies_field(nnx,nnode,x,y,1)
            ina1 = indicies_field(nnx,nnode,x,y,2)
            ina2 = indicies_field(nnx,nnode,x,y,3)
            ina3 = indicies_field(nnx,nnode,x,y,4)
            ina4 = indicies_field(nnx,nnode,x,y,5)
            ina5 = indicies_field(nnx,nnode,x,y,6)
            ina6 = indicies_field(nnx,nnode,x,y,7)
            ina7 = indicies_field(nnx,nnode,x,y,8)
            ina8 = indicies_field(nnx,nnode,x,y,9)

            tl0 = t2[ina0]
            tl1 = t2[ina1]
            tl2 = t2[ina2]
            tl3 = t2[ina3]
            tl4 = t2[ina4]
            tl5 = t2[ina5]
            tl6 = t2[ina6]
            tl7 = t2[ina7]
            tl8 = t2[ina8]

            ct = tl0 + tl1 + tl2 + tl3 + tl4 + tl5 + tl6 + tl7 + tl8

            temp[in] = ct

            frilc = 0.0
            frjlc = betah * (ct - tmid) * Fo
            fri[in] = frilc
            frj[in] = frjlc

            fl0 = f2[ina0]
            fl1 = f2[ina1]
            fl2 = f2[ina2]
            fl3 = f2[ina3]
            fl4 = f2[ina4]
            fl5 = f2[ina5]
            fl6 = f2[ina6]
            fl7 = f2[ina7]
            fl8 = f2[ina8]

            cf = fl0 + fl1 + fl2 + fl3 + fl4 + fl5 + fl6 + fl7 + fl8
            cui = (fl1 + fl5 + fl8 - fl3 - fl6 - fl7)
            cuj = (fl2 + fl5 + fl6 - fl4 - fl7 - fl8)

            rho[in] = cf
            invrho = 1.0 / (cf)

            ui[in] = invrho * (cui + 0.5 * frilc) * dt
            uj[in] = invrho * (cuj + 0.5 * frjlc) * dt 

            pres[in] = 1.0 + 0.5 * (frilc * ui[in] + frjlc * uj[in])

        end 
    end

    xm=nnx-1

    for y = 1 : nny

        inb = indicies_scalar(nnx,1,y)
        inbp = indicies_scalar(nnx,2,y)
        int = indicies_scalar(nnx,nnx,y)
        intm = indicies_scalar(nnx,xm,y)

        ui[inb]=uwxb
        uj[inb]=uwyb
        tempp = temp[inbp]
        temp[inb] = tempp

        ui[int]=uwxt
        uj[int]=uwyt
        tempm = temp[intm]
        temp[int]= tempm
    end

    for x = 1 : nnx

        inl = indicies_scalar(nnx,x,1)
        inr = indicies_scalar(nnx,x,nny)

        ui[inl]=uwxl
        uj[inl]=uwyl
        temp[inl]=thot
        
        ui[inr]=uwxr
        uj[inr]=uwyr
        temp[inr]=tcold
    end
end
