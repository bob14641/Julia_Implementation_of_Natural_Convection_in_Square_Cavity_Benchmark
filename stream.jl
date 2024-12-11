function stream(nnx,nny,nnode,ndir,vpi,solidnode,solidstr,f1,f2,t1,t2,temp,thot,tcold,phi,phip,phim,invphip,omegat,uwxr,uwxl,uwxb,uwxt,uwyr,uwyl,uwyb,uwyt,cs2,dt,fbx,fby)

for x = 1 : nnx
    for y = 1 : nny 
        for d = 1 : ndir
            ys = y - vpi[d,1]
            xs = x - vpi[d,2]
            in = indicies_scalar(nnx,x,y)
            inp = indicies_scalar(nnx,xs,ys)
            ina = indicies_field(nnx,nnode,x,y,d)
            inap = indicies_field(nnx,nnode,xs,ys,d)
            if (xs >= 1) && (xs <= nnx)
                if  (ys >= 1) && (ys <= nny)
                    if (solidstr[in] == 0) && (solidstr[inp] == 0) 
                    f2[ina] = f1[inap]
                    t2[ina] = t1[inap]
                    end
                end
            end
        end
    end
end

for x = 1:nnx
    for y = 1 : nny
        in = indicies_scalar(nnx,x,y)
        in0 = indicies_field(nnx,nnode,x,y,1)
        in1 = indicies_field(nnx,nnode,x,y,2)
        in2 = indicies_field(nnx,nnode,x,y,3)
        in3 = indicies_field(nnx,nnode,x,y,4)
        in4 = indicies_field(nnx,nnode,x,y,5)
        in5 = indicies_field(nnx,nnode,x,y,6)
        in6 = indicies_field(nnx,nnode,x,y,7)
        in7 = indicies_field(nnx,nnode,x,y,8)
        in8 = indicies_field(nnx,nnode,x,y,9)

        if (solidnode[in] == 1) #Bottom Wall

            #===================================================================================#
            #Flow Part
            #===================================================================================#

            b = 1.0 / (1.0 - uwyb) 

            fs0 = f2[in0]
            fs1 = f2[in1]
            fs3 = f2[in3]
            fs4 = f2[in4]
            fs7 = f2[in7]
            fs8 = f2[in8]

            rwb = b * (fs0 + fs1 + fs3 + 2.0 * (fs4 + fs7 + fs8) - 0.5 * fby * dt)

            C10b = rwb * (uwxb) - 0.5 * fbx * dt - (fs1 - fs3 - fs7 + fs8) 
            C01b = rwb * (uwyb) - 0.5 * fby * dt - (-fs4 - fs7 -fs8) 
            C2002b = (2.0 * phi * cs2 * rwb) + (phip * rwb * uwxb * uwxb) + (phim * rwb * uwyb * uwyb) - (phip * fbx * uwxb * dt) - (phim * fby * uwyb * dt) - (phi * (fs1 + fs3 + fs4 + 2.0 * fs7 + 2.0 * fs8) + (fs1 + fs3 - fs4))

            fs2 = invphip * (2.0 * phi * C01b - C2002b)
            fs5 = 0.5 * invphip * (C01b + C10b + C2002b + (-C01b + C10b) * phi) 
            fs6 = -0.5 * invphip * (C10b - C01b - C2002b + (C01b + C10b) * phi) 

            f2[in2] = fs2
            f2[in5] = fs5 
            f2[in6] = fs6

            #===================================================================================#
            #Thermal Part
            #===================================================================================#

            my = y - 1
            py = 1 + y

            inp0 = indicies_field(nnx,nnode,x,py,1)
            inp1 = indicies_field(nnx,nnode,x,py,2)
            inp2 = indicies_field(nnx,nnode,x,py,3)
            inp3 = indicies_field(nnx,nnode,x,py,4)
            inp4 = indicies_field(nnx,nnode,x,py,5)
            inp5 = indicies_field(nnx,nnode,x,py,6)
            inp6 = indicies_field(nnx,nnode,x,py,7)
            inp7 = indicies_field(nnx,nnode,x,py,8)
            inp8 = indicies_field(nnx,nnode,x,py,9)
    
            inm0 = indicies_field(nnx,nnode,x,my,1)
            inm1 = indicies_field(nnx,nnode,x,my,2)
            inm2 = indicies_field(nnx,nnode,x,my,3)
            inm3 = indicies_field(nnx,nnode,x,my,4)
            inm4 = indicies_field(nnx,nnode,x,my,5)
            inm5 = indicies_field(nnx,nnode,x,my,6)
            inm6 = indicies_field(nnx,nnode,x,my,7)
            inm7 = indicies_field(nnx,nnode,x,my,8)
            inm8 = indicies_field(nnx,nnode,x,my,9)

            ts0 = t2[in0]
            ts1 = t2[in1]
            ts3 = t2[in3]
            ts4 = t2[in4]
            ts7 = t2[in7]
            ts8 = t2[in8]

            tsp0 = t2[inp0]
            tsp1 = t2[inp1]
            tsp3 = t2[inp3]
            tsp4 = t2[inp4]
            tsp7 = t2[inp7]
            tsp8 = t2[inp8]

            tsm0 = t2[inm0]
            tsm1 = t2[inm1]
            tsm3 = t2[inm3]
            tsm4 = t2[inm4]
            tsm7 = t2[inm7]
            tsm8 = t2[inm8]

            tlcb = b * (ts0 + ts1 + ts3 + 2.0 * (ts4 + ts7 + ts8))
            tlcbp = b * (tsp0 + tsp1 + tsp3 + 2.0 * (tsp4 + tsp7 + tsp8))
            tlcbm = b * (tsm0 + tsm1 + tsm3 + 2.0 * (tsm4 + tsm7 + tsm8))

            dtx = 0.5 * (tlcbp - tlcbm)

            E00b = tlcb - (ts0 + ts1 + ts3 + ts4 + ts7 + ts8)
            E10b = -(cs2 * dtx / omegat) - (ts1 - ts3 - ts7 + ts8)
            E20b = cs2 * tlcb + tlcb * (uwxb * uwxb) - (ts1 + ts3 + ts7 + ts8)

            ts2 = E00b - E20b
            ts5 = 0.5 * (E10b + E20b)
            ts6 = 0.5 * (E20b - E10b)

            t2[in2] = ts2
            t2[in5] = ts5
            t2[in6] = ts6

        elseif (solidnode[in] == 2) #Left Wall

            #===================================================================================#
            #Flow Part
            #===================================================================================#

            l = 1.0 / (1.0 - uwxl) 

            fs0 = f2[in0]
            fs2 = f2[in2]
            fs3 = f2[in3]
            fs4 = f2[in4]
            fs6 = f2[in6]
            fs7 = f2[in7]

            rwl = l * (fs0 + fs2 + fs4 + 2.0 * (fs3 + fs6 + fs7) - 0.5 * fbx * dt)

            C10l = rwl * (uwxl) - (-fs3 - fs6 - fs7) - 0.5 * fbx * dt
            C01l = rwl * (uwyl) - (fs2 - fs4 + fs6 - fs7) - 0.5 * fby * dt
            C2002l = (2.0 * phi * cs2 * rwl) + (phim * rwl * uwxl * uwxl) + (phip * rwl * uwyl * uwyl) - (phim * fbx * uwxl * dt) - (phip * fby * uwyl * dt) - (phi * (fs2 + fs3 + fs4 + 2.0 * fs6 + 2.0 * fs7) - (-fs2 + fs3 - fs4))

            fs1 =  invphip * (2.0 * phi * C10l - C2002l)
            fs5 = 0.5 * invphip * (C01l + C10l + C2002l + (C01l - C10l) * phi)
            fs8 = - 0.5 * invphip * (C01l - C10l - C2002l + (C01l + C10l) * phi)

            f2[in1] = fs1
            f2[in5] = fs5 
            f2[in8] = fs8

            #===================================================================================#
            #Thermal Part
            #===================================================================================#

            tlcl = thot
            
            ts0 = t2[in0]
            ts2 = t2[in2]
            ts3 = t2[in3]
            ts4 = t2[in4]
            ts6 = t2[in6]
            ts7 = t2[in7]

            E00l = tlcl - (ts0 + ts2 + ts3 + ts4 + ts6 + ts7)
            E01l = tlcl * uwyl - (ts2 - ts4 + ts6 - ts7)
            E02l = cs2 * tlcl + tlcl * (uwyl * uwyl) - (ts2 + ts4 + ts6 + ts7)

            ts1 = E00l - E02l
            ts5 = 0.5 * (E01l + E02l)
            ts8 = 0.5 * (E02l - E01l)

            t2[in1] = ts1
            t2[in5] = ts5
            t2[in8] = ts8

        elseif (solidnode[in] == 3) #Top Wall

            #===================================================================================#
            #Flow Part
            #===================================================================================#

            t = 1.0 / (1.0 + uwyt) 

            fs0 = f2[in0]
            fs1 = f2[in1]
            fs2 = f2[in2]
            fs3 = f2[in3]
            fs5 = f2[in5]
            fs6 = f2[in6]
            
            rwt = t * (fs0 + fs1 + fs3 + 2.0 * (fs2 + fs5 + fs6) + 0.5 * fby * dt)

            C10t = rwt * (uwxt) - (fs1 - fs3 + fs5 - fs6) - 0.5 * fbx * dt
            C01t = rwt * (uwyt) - (fs2 + fs5 + fs6) - 0.5 * fby * dt
            C2002t = (2.0 * phi * cs2 * rwt) + (phip * rwt * uwxt * uwxt) + (phim * rwt * uwyt * uwyt) - (phip * fbx * uwxt * dt) - (phim * fby * uwyt * dt) - (phi * (fs1 + fs2 + fs3 + 2.0 * fs5 + 2.0 * fs6) + (fs1 - fs2 + fs3))

            fs4 = invphip * (-2.0 * phi * C01t - C2002t)
            fs7 = -0.5 * invphip * (C01t + C10t - C2002t + (-C01t + C10t) * phi)
            fs8 = 0.5 * invphip * (C10t - C01t + C2002t + (C01t + C10t) * phi)

            f2[in4] = fs4
            f2[in7] = fs7 
            f2[in8] = fs8

            #===================================================================================#
            #Thermal Part
            #===================================================================================#

            my = y - 1
            py = 1 + y

            inp0 = indicies_field(nnx,nnode,x,py,1)
            inp1 = indicies_field(nnx,nnode,x,py,2)
            inp2 = indicies_field(nnx,nnode,x,py,3)
            inp3 = indicies_field(nnx,nnode,x,py,4)
            inp4 = indicies_field(nnx,nnode,x,py,5)
            inp5 = indicies_field(nnx,nnode,x,py,6)
            inp6 = indicies_field(nnx,nnode,x,py,7)
            inp7 = indicies_field(nnx,nnode,x,py,8)
            inp8 = indicies_field(nnx,nnode,x,py,9)
    
            inm0 = indicies_field(nnx,nnode,x,my,1)
            inm1 = indicies_field(nnx,nnode,x,my,2)
            inm2 = indicies_field(nnx,nnode,x,my,3)
            inm3 = indicies_field(nnx,nnode,x,my,4)
            inm4 = indicies_field(nnx,nnode,x,my,5)
            inm5 = indicies_field(nnx,nnode,x,my,6)
            inm6 = indicies_field(nnx,nnode,x,my,7)
            inm7 = indicies_field(nnx,nnode,x,my,8)
            inm8 = indicies_field(nnx,nnode,x,my,9)

            ts0 = t2[in0]
            ts1 = t2[in1]
            ts2 = t2[in2]
            ts3 = t2[in3]
            ts5 = t2[in5]
            ts6 = t2[in6]
            
            tsp0 = t2[inp0]
            tsp1 = t2[inp1]
            tsp2 = t2[inp2]
            tsp3 = t2[inp3]
            tsp5 = t2[inp5]
            tsp6 = t2[inp6]

            tsm0 = t2[inm0]
            tsm1 = t2[inm1]
            tsm2 = t2[inm2]
            tsm3 = t2[inm3]
            tsm5 = t2[inm5]
            tsm6 = t2[inm6]

            tlct = t * (ts0 + ts1 + ts3 + 2.0 * (ts2 + ts5 + ts6))
            tlctp = t * (tsp0 + tsp1 + tsp3 + 2.0 * (tsp2 + tsp5 + tsp6))
            tlctm = t * (tsm0 + tsm1 + tsm3 + 2.0 * (tsm2 + tsm5 + tsm6))

            dtx = 0.5 * (tlctp - tlctm)

            E00t = tlct - (ts0 + ts1 + ts2 + ts3 + ts5 + ts6)
            E10t = -(cs2 * dtx / omegat) - (ts1 - ts3 + ts5 - ts6)
            E20t = cs2 * tlct + tlct * (uwxt * uwxt) - (ts1 + ts3 + ts5 + ts6)

            ts4 = E00t - E20t
            ts7 = 0.5 * (E20t - E10t)
            ts8 = 0.5 * (E20t + E10t)

            t2[in4] = ts4
            t2[in7] = ts7
            t2[in8] = ts8

        elseif (solidnode[in] == 4) #Right Wall

            #===================================================================================#
            #Flow Part
            #===================================================================================#

            r = 1.0 / (1.0 + uwxr) 

            fs0 = f2[in0]
            fs2 = f2[in2]   
            fs4 = f2[in4]

            fs1 = f2[in1]
            fs5 = f2[in5]
            fs8 = f2[in8]

            rwr = r * (fs0 + fs2 + fs4 + 2.0 * (fs1 + fs5 + fs8) + 0.5 * fbx * dt)

            C10r = rwr * (uwxr) - (fs1 + fs5 + fs8) - 0.5 * fbx * dt
            C01r = rwr * (uwyr) - (fs2 - fs4 + fs5 - fs8) - 0.5 * fby * dt
            C2002r = (2.0 * phi * cs2 * rwr) + (phim * rwr * uwxr * uwxr) + (phip * rwr * uwyr * uwyr) - (phim * fbx * uwxr * dt) - (phip * fby * uwyr * dt) - (phi*(fs1 + fs2 + fs4 + 2.0*fs5 + 2.0*fs8) - (fs1 - fs2 - fs4))

            fs3 = invphip * ((-2.0 * phi * C10r) - C2002r)
            fs6 = 0.5 * invphip * (C01r - C10r + C2002r + (C01r + C10r) * phi)
            fs7 = -0.5 * invphip * (C01r + C10r - C2002r + (C01r - C10r) * phi)

            f2[in3] = fs3
            f2[in6] = fs6
            f2[in7] = fs7

            #===================================================================================#
            #Thermal Part
            #===================================================================================#

            tlcr = tcold

            ts0 = t2[in0]
            ts1 = t2[in1]
            ts2 = t2[in2]
            ts4 = t2[in4]
            ts5 = t2[in5]
            ts8 = t2[in8]

            E00r = tlcr - (ts0 + ts1 + ts2 + ts4 + ts5 + ts8)
            E01r = tlcr * uwyr - (ts2 - ts4 + ts5 - ts8)
            E02r = cs2 * tlcr + tlcr * (uwyr * uwyr) - (ts2 + ts4 + ts5 + ts8)

            ts3 = E00r - E02r
            ts6 = 0.5 * (E01r + E02r)
            ts7 = 0.5 * (E02r - E01r)

            t2[in3] = ts3
            t2[in6] = ts6
            t2[in7] = ts7

        elseif (solidnode[in] == 5) #Bottom-Left Corner

            #===================================================================================#
            #Flow Part
            #===================================================================================#

            fs0 = f2[in0]
            fs3 = f2[in3]
            fs4 = f2[in4]
            fs7 = f2[in7]

            rw7 = - (fs0 + 2.0 * fs3 + 2.0 * fs4 + 4.0 * fs7) / ((uwxb) + (uwyb) - (uwxb * uwyb) - 1.0)

            C10bl = fs3 + fs7 + rw7 * (uwxb) 
            C01bl = fs4 + fs7 + rw7 * (uwyb)
            C2sbl = (2.0 * cs2 * rw7) - fs4 - 2.0 * fs7 - fs3 + rw7 * ((uwxb * uwxb) + (uwyb * uwyb))
            C2dbl = fs4 - fs3 + rw7 * ((uwxb * uwxb) + (uwyb * uwyb))
            C11bl = (rw7 * (uwxb * uwxb * uwyb * uwyb)) - fs7
    
            fs1 = C01bl + C10bl - C11bl + 0.5 * (C2dbl - C2sbl)
            fs2 = C01bl + C10bl - C11bl + 0.5 * (-C2dbl - C2sbl)
            fs5 = C11bl + 0.5 * (-C10bl - C01bl + C2sbl)
            fs6 = 0.5 * (-C10bl) + 0.25 * (C2dbl + C2sbl)
            fs8 = 0.5 * (-C01bl) + 0.25 * (-C2dbl + C2sbl)

            f2[in1] = fs1 
            f2[in2] = fs2
            f2[in5] = fs5
            f2[in6] = fs6
            f2[in8] = fs8 

            #===================================================================================#
            #Thermal Part
            #===================================================================================#

            ts0=t2[in0]
            ts3=t2[in3]
            ts4=t2[in4]
            ts7=t2[in7]

            twbl = thot 

            gx10 = (ts0 + 2.0 * ts3 + 2.0 * ts4 + 4.0 * ts7 - twbl) / (cs2 / omegat)
            gy01 = 0.0
            gx20 = 0.0
            gy02 = 0.0
            gx11 = 0.0
            gy11 = 0.0

            E10bl = twbl * uwxb - (cs2 * gx10 / omegat) - (-ts3 - ts7)
            E01bl = twbl * uwyb - gy01 - (-ts4 - ts7)
            E20bl = twbl * cs2 + twbl * (uwxb * uwxb) - gx20 - (ts3 + ts7)
            E02bl = twbl * cs2 + twbl * (uwyb * uwyb) - gy02 - (ts4 + ts7)
            E11bl = twbl * uwxb * uwyb - gx11 - gy11 - ts7 

            ts1 = E10bl + E01bl - E02bl - E11bl 
            ts2 = E10bl + E01bl - E20bl - E11bl 
            ts5 = E11bl + 0.5 * (E20bl + E02bl - E10bl - E01bl)
            ts6 = 0.5 * (E20bl - E10bl)
            ts8 = 0.5 * (E02bl - E01bl)

            t2[in1] = ts1 
            t2[in2] = ts2
            t2[in5] = ts5
            t2[in6] = ts6
            t2[in8] = ts8 
 
        elseif (solidnode[in] == 6) #Bottom-Right Corner

            #===================================================================================#
            #Flow Part
            #===================================================================================#

            fs0 = f2[in0]
            fs1 = f2[in1]
            fs4 = f2[in4]
            fs8 = f2[in8]

            rw8 = (fs0 + 2.0 * fs1 + 2.0 * fs4 + 4.0 * fs8) / ((uwxb) - (uwyb) - (uwxb * uwyb) + 1.0)

            C10br = -fs1 - fs8 + rw8 * (uwxb) 
            C01br = fs4 + fs8 + rw8 * (uwyb)
            C2sbr = (2.0 * cs2 * rw8) - fs4 - 2.0 * fs8 - fs1 + rw8 * ((uwxb * uwxb) + (uwyb  * uwyb))
            C2dbr = fs4 - fs1 + rw8 * ((uwxb * uwxb) + (uwyb * uwyb))
            C11br = (rw8 * (uwxb * uwxb) * (uwyb * uwyb)) + fs8
    
            fs2 = C01br - C10br + C11br + 0.5 * (-C2dbr - C2sbr)
            fs3 = C01br - C10br + C11br + 0.5 * (C2dbr - C2sbr)
            fs5 = 0.5 * (C10br) + 0.25 * (C2dbr + C2sbr)
            fs6 = -C11br + 0.5 * (C10br - C01br + C2sbr)
            fs7 = 0.5 * (-C01br) + 0.25 * (-C2dbr + C2sbr)

            f2[in2] = fs2 
            f2[in3] = fs3
            f2[in5] = fs5
            f2[in6] = fs6
            f2[in7] = fs7

            #===================================================================================#
            #Thermal Part
            #===================================================================================#

            ts0=t2[in0]
            ts1=t2[in1]
            ts4=t2[in4]
            ts8=t2[in8]
            twbr = tcold

            gx10 = (twbr - (ts0 + 2.0 * ts1 + 2.0 * ts4 + 4.0 * ts8)) / (cs2 / omegat)
            gy01 = 0.0
            gx20 = 0.0
            gy02 = 0.0
            gx11 = 0.0
            gy11 = 0.0

            E10br = twbr * uwxb -(cs2 * gx10 / omegat) - (ts1 + ts8)
            E01br = twbr * uwyb - gy01 - (-ts4 - ts8)
            E20br = twbr * cs2 + twbr * (uwxb * uwxb) - gx20 - (ts1 + ts8)
            E02br = twbr * cs2 + twbr * (uwyb * uwyb) - gy02 - (ts4 + ts8)
            E11br = twbr * uwxb * uwyb - gx11 - gy11 + ts8 

            ts2 = E01br + E11br - E10br - E20br 
            ts3 = E01br + E11br - E10br - E02br 
            ts5 =  0.5 * (E20br + E10br)
            ts6 = -E11br + 0.5 * (E20br + E02br + E10br - E01br)
            ts7 = 0.5 * (E02br - E01br)

            t2[in2] = ts2 
            t2[in3] = ts3
            t2[in5] = ts5
            t2[in6] = ts6
            t2[in7] = ts7

        elseif (solidnode[in] == 7) #Top-Left Corner

            #===================================================================================#
            #Flow Part
            #===================================================================================#

            fs0 = f2[in0]
            fs2 = f2[in2]
            fs3 = f2[in3]
            fs6 = f2[in6]

            rw5 = - (fs0 + 2.0 * fs2 + 2.0 * fs3 + 4.0 * fs6) / ((uwxt) - (uwyt) + (uwxt * uwyt) - 1.0)

            C10tl = fs3 + fs6 + rw5 * (uwxt) 
            C01tl = -fs2 - fs6 + rw5 * (uwyt)
            C2stl = (2.0 * cs2 * rw5) - fs3 - 2.0 * fs6 - fs2 + rw5 * ((uwxt * uwxt) + (uwyt * uwyt))
            C2dtl = fs2 - fs3 + rw5 * ((uwxt * uwxt) - (uwyt * uwyt))
            C11tl = (rw5 * (uwxt * uwxt) * (uwyt * uwyt)) + fs6
    
            fs1 = C10tl - C01tl + C11tl + 0.5 * (C2dtl - C2stl)
            fs4 = C10tl - C01tl + C11tl + 0.5 * (-C2dtl - C2stl)
            fs5 = 0.5 * (C01tl) + 0.25 * (-C2dtl + C2stl) 
            fs7 = 0.5 * (-C10tl) + 0.25 * (C2dtl + C2stl)
            fs8 = -C11tl + 0.5 * (C01tl - C10tl + C2stl)

            f2[in1] = fs1
            f2[in4] = fs4
            f2[in5] = fs5 
            f2[in7] = fs7
            f2[in8] = fs8

            #===================================================================================#
            #Thermal Part
            #===================================================================================#

            ts0=t2[in0]
            ts2=t2[in2]
            ts3=t2[in3]
            ts6=t2[in6]
            
            twtl = thot 

            gx10 = (ts0 + 2.0 * ts2 + 2.0 * ts3 + 4.0 * ts6 - twtl) / (cs2 / omegat)
            gy01 = 0.0
            gx20 = 0.0
            gy02 = 0.0
            gx11 = 0.0
            gy11 = 0.0

            E10tl = twtl * uwxt - (cs2 * gx10 / omegat) - (-ts3 - ts6)
            E01tl = twtl * uwyt - gy01 - (ts2 + ts6)
            E20tl = twtl * cs2 + twtl * (uwxt * uwxt) - gx20 - (ts3 + ts6)
            E02tl = twtl * cs2 + twtl * (uwyt * uwyt) - gy02 - (ts2 + ts6)
            E11tl = twtl * uwxt * uwyt - gx11 - gy11 + ts6 

            ts1 = E10tl - E01tl - E02tl + E11tl 
            ts4 = E10tl - E01tl - E20tl + E11tl 
            ts5 = 0.5 * (E01tl + E02tl)
            ts7 = 0.5 * (E20tl - E10tl)
            ts8 = -E11tl + 0.5 * (E20tl + E02tl - E10tl + E01tl)
            
            t2[in1] = ts1
            t2[in4] = ts4
            t2[in5] = ts5 
            t2[in7] = ts7
            t2[in8] = ts8

        elseif (solidnode[in] == 8) #Top-Right Corner

            #===================================================================================#
            #Flow Part
            #===================================================================================#

            fs0 = f2[in0]
            fs1 = f2[in1]
            fs2 = f2[in2]
            fs5 = f2[in5]

            rw6 = (fs0 + 2.0 * fs1 + 2.0 * fs2 + 4.0 * fs5) / ((uwxt) + (uwyt) + (uwxt * uwyt) + 1.0)

            C10tr = - fs5 - fs1 + rw6 * (uwxt) 
            C01tr = - fs5 - fs2 + rw6 * (uwyt)
            C2str = (2.0 * cs2 * rw6) - fs2 - 2.0 * fs5 - fs1 + rw6 * ((uwxt * uwxt) + (0uwyt * uwyt))
            C2dtr = fs2 - fs1 + rw6 * ((uwxt * uwxt) - (uwyt * uwyt))
            C11tr = (rw6 * (uwxt * uwxt) * (uwyt * uwyt)) - fs5
    
            fs3 = -C10tr - C01tr - C11tr + 0.5 * (C2dtr - C2str)
            fs4 = -C10tr - C01tr - C11tr + 0.5 * (- C2dtr - C2str)
            fs6 = 0.5 * (C01tr) + 0.25 * (-C2dtr + C2str) 
            fs7 = C11tr + 0.5 * (C01tr + C10tr + C2str)
            fs8 = 0.5 * (C10tr) + 0.25 * (C2dtr + C2str)

            f2[in3] = fs3 
            f2[in4] = fs4
            f2[in6] = fs6
            f2[in7] = fs7
            f2[in8] = fs8    
            
            #===================================================================================#
            #Thermal Part
            #===================================================================================#

            ts0 = t2[in0]
            ts1 = t2[in1]
            ts2 = t2[in2]
            ts5 = t2[in5]

            twtr = tcold

            gx10 = (twtr - (ts0 + 2.0 * ts1 + 2.0 * ts2 + 4.0 * ts5)) / (cs2 / omegat)
            gy01 = 0.0
            gx20 = 0.0
            gy02 = 0.0
            gx11 = 0.0
            gy11 = 0.0

            E10tr = twtr * uwxt - (cs2 * gx10 / omegat) - (ts1 + ts5)
            E01tr = twtr * uwyt - gy01 - (ts2 + ts5)
            E20tr = twtr * cs2 + twtr * (uwxt * uwxt) - gx20 - (ts1 + ts5)
            E02tr = twtr * cs2 + twtr * (uwyt * uwyt) - gy02 - (ts2 + ts5)
            E11tr = twtr * uwxt * uwyt - gx11 - gy11 - ts5

            ts3 = -E10tr - E01tr - E02tr - E11tr 
            ts4 = -E10tr - E01tr - E20tr - E11tr 
            ts6 = 0.5 * (E01tr + E02tr)
            ts7 = E11tr + 0.5 * (E20tr + E02tr + E10tr + E01tr)
            ts8 = 0.5 * (E20tr + E10tr)

            t2[in3] = ts3 
            t2[in4] = ts4
            t2[in6] = ts6
            t2[in7] = ts7
            t2[in8] = ts8
        end
    end 
end

return f1

end