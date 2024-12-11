function collide_CM_thermal(tauf,taut,nnx,nny,nnode,solidnode,ui,uj,rho,f2,t2,cs2,temp,pres,fri,frj,tmid,Fo,betah)

wf0 = 1.0
wf1 = 1.0
wf2 = 1.0
wf3 = 1.0 
wf4 = 1.0 / tauf
wf5 = 1.0 / tauf
wf6 = 1.0
wf7 = 1.0
wf8 = 1.0

wt0 = 1.0
wt1 = 1.0 / taut 
wt2 = 1.0 / taut
wt3 = 1.0
wt4 = 1.0
wt5 = 1.0 
wt6 = 1.0
wt7 = 1.0 
wt8 = 1.0
    
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

        templc = temp[in]
        preslc = pres[in]
        uilc = ui[in]
        ujlc = uj[in]
        uilc2 = uilc * uilc
        ujlc2 = ujlc * ujlc
        uiujlc = uilc * ujlc
        rholc = rho[in]
        invrholc = 1.0 / rholc
        uilctwo = 2.0 * uilc
        ujlctwo = 2.0 * ujlc
        uiujlctwo = 2.0 * uiujlc
        frilc = fri[in]
        frjlc = frj[in]

##---------------------------------------------------------------------------------------------------------   
## Fluid Motion Part of Collision
##--------------------------------------------------------------------------------------------------------- 
            
        fl0 = f2[ina0]
        fl1 = f2[ina1]
        fl2 = f2[ina2]
        fl3 = f2[ina3]
        fl4 = f2[ina4]
        fl5 = f2[ina5]
        fl6 = f2[ina6]
        fl7 = f2[ina7]
        fl8 = f2[ina8]

        f13 = fl1 + fl3
        f24 = fl2 + fl4
        f58a = fl5 + fl6 + fl7 + fl8
        f58b = fl5 - fl6 - fl7 + fl8
        f58c = fl5 + fl6 - fl7 - fl8

        Kp00 = fl0 + f13 + f24 + f58a
        Kp10 = fl1 - fl3 + f58b
        Kp01 = fl2 - fl4 + f58c 
        Kp20 = f13 + f58a
        Kp02 = f24 + f58a
        Kp11 = fl5 - fl6 + fl7 - fl8
        Kp21 = f58c
        Kp12 = f58b
        Kp22 = f58a
            
        K00 = Kp00
        K10 = Kp10 - uilc * Kp00
        K01 = Kp01 - ujlc * Kp00
        K20 = Kp20 - uilctwo * Kp10 + uilc2 * Kp00
        K02 = Kp02 - ujlctwo * Kp01 + ujlc2 * Kp00
        K11 = Kp11 - ujlc * Kp10 - uilc * Kp01 + uiujlc * Kp00 
        K21 = Kp21 - uilctwo * Kp11 + uilc2 * Kp01 - ujlc * Kp20 + uiujlctwo * Kp10 - uilc2 * ujlc * Kp00
        K12 = Kp12 - ujlctwo * Kp11 + ujlc2 * Kp10 - uilc * Kp02 + uiujlctwo * Kp01 - uilc * ujlc2 * Kp00
        K22 = Kp22 - uilctwo * Kp12 + uilc2 * Kp02 - ujlctwo * Kp21 + 4.0 * uiujlc * Kp11 - uilc2 * ujlctwo * Kp01 + ujlc2 * Kp20 - uilctwo * ujlc2 * Kp10 + uilc2 * ujlc2 * Kp00
    
        K2s = K20 + K02
        K2d = K20 - K02

        Keq_00 = rholc
        Keq_10 = 0.0
        Keq_01 = 0.0
        Keq_20 = cs2 * rholc
        Keq_02 = cs2 * rholc
        Keq_11 = 0.0
        Keq_21 = 0.0
        Keq_12 = 0.0
        Keq_22 = cs2 * cs2 * rholc

        Keq_2s = Keq_20 + Keq_02
        Keq_2d = Keq_20 - Keq_02

        S00 = 0.0 
        S10 = frilc 
        S01 = frjlc
        S20 = 2.0 * frilc * uilc 
        S02 = 2.0 * frjlc * ujlc
        S11 = frilc * ujlc + frjlc * uilc
        S21 = 0.0
        S12 = 0.0

        S2s = S20 + S02
        S2d = S20 - S02
    
        K00 = K00 + (wf0 * (Keq_00 - K00) + (1.0 - 0.5 * wf0) * S00)
        K10 = K10 + (wf1 * (Keq_10 - K10) + (1.0 - 0.5 * wf1) * S10)
        K01 = K01 + (wf2 * (Keq_01 - K01) + (1.0 - 0.5 * wf2) * S01)
        K2s = K2s + (wf3 * (Keq_2s - K2s) + (1.0 - 0.5 * wf3) * S2s)
        K2d = K2d + (wf4 * (Keq_2d - K2d) + (1.0 - 0.5 * wf4) * S2d)
        K11 = K11 + (wf5 * (Keq_11 - K11) + (1.0 - 0.5 * wf5) * S11)
        K21 = K21 + wf6 * (Keq_21 - K21)
        K12 = K12 + wf7 * (Keq_12 - K12)

        K20 = 0.5 * (K2s + K2d)
        K02 = 0.5 * (K2s - K2d)

        Keq_22 = (invrholc) * (K20 * K02 + 2.0 * K11 * K11)
        K22 = K22 + wf8 * (Keq_22 - K22)
    
        Kp00 = K00
        Kp10 = K10 + uilc * K00
        Kp01 = K01 + ujlc * K00
        Kp20 = K20 + uilctwo * K10 + uilc2 * K00
        Kp02 = K02 + ujlctwo * K01 + ujlc2 * K00
        Kp11 = K11 + uilc * K01 + ujlc * K10 + uilc * ujlc * K00
        Kp21 = K21 + uilctwo * K11 + ujlc * K20 + uilc2 * K01 + uilctwo * ujlc * K10 + uilc2 * ujlc * K00
        Kp12 = K12 + ujlctwo * K11 + uilc * K02 + uilctwo * ujlc * K01 + ujlc2 * K10 + uilc * ujlc2 * K00
        Kp22 = K22 + uilctwo * K12 + ujlctwo * K21 + 4.0 * uiujlc * K11 + ujlc2 * K20 + uilc2 * K02 + uilc2 * ujlctwo * K01 + uilctwo * ujlc2 * K10 + uilc2 * ujlc2 * K00
    
        f2[ina0] = Kp00 - Kp20 - Kp02 + Kp22
        f2[ina1] = 0.5 * (Kp10 + Kp20 - Kp12 - Kp22)
        f2[ina2] = 0.5 * (Kp01 + Kp02 - Kp21 - Kp22)
        f2[ina3] = 0.5 * (-Kp10 + Kp20 + Kp12 - Kp22)
        f2[ina4] = 0.5 * (-Kp01 + Kp02 + Kp21 - Kp22)
        f2[ina5] = 0.25 * (Kp11 + Kp21 + Kp12 + Kp22)
        f2[ina6] = 0.25 * (-Kp11 + Kp21 - Kp12 + Kp22)
        f2[ina7] = 0.25 * (Kp11 - Kp21 - Kp12 + Kp22)
        f2[ina8] = 0.25 * (-Kp11 - Kp21 + Kp12 + Kp22)

##---------------------------------------------------------------------------------------------------------   
## Thermal Part of Collision
##--------------------------------------------------------------------------------------------------------- 
                
        tl0 = t2[ina0]
        tl1 = t2[ina1]
        tl2 = t2[ina2]
        tl3 = t2[ina3]
        tl4 = t2[ina4]
        tl5 = t2[ina5]
        tl6 = t2[ina6]
        tl7 = t2[ina7]
        tl8 = t2[ina8]

        t13 = tl1 + tl3
        t24 = tl2 + tl4
        t58a = tl5 + tl6 + tl7 + tl8
        t58b = tl5 - tl6 - tl7 + tl8
        t58c = tl5 + tl6 - tl7 - tl8

        Xp00 = tl0 + t13 + t24 + t58a
        Xp10 = tl1 - tl3 + t58b
        Xp01 = tl2 - tl4 + t58c 
        Xp20 = t13 + t58a
        Xp02 = t24 + t58a
        Xp11 = tl5 - tl6 + tl7 - tl8
        Xp21 = t58c
        Xp12 = t58b
        Xp22 = t58a 
                
        X00 = Xp00
        X10 = Xp10 - uilc * Xp00
        X01 = Xp01 - ujlc * Xp00
        X20 = Xp20 - uilctwo * Xp10 + uilc2 * Xp00
        X02 = Xp02 - ujlctwo * Xp01 + ujlc2 * Xp00
        X11 = Xp11 - ujlc * Xp10 - uilc * Xp01 + uiujlc * Xp00 
        X21 = Xp21 - uilctwo * Xp11 + uilc2 * Xp01 - ujlc * Xp20 + uiujlctwo * Xp10 - uilc2 * ujlc * Xp00
        X12 = Xp12 - ujlctwo * Xp11 + ujlc2 * Xp10 - uilc * Xp02 + uiujlctwo * Xp01 - uilc * ujlc2 * Xp00
        X22 = Xp22 - uilctwo * Xp12 + uilc2 * Xp02 - ujlctwo * Xp21 + 4.0 * uiujlc * Xp11 - uilc2 * ujlctwo * Xp01 + ujlc2 * Xp20 - uilctwo * ujlc2 * Xp10 + uilc2 * ujlc2 * Xp00
        
        X2s = X20 + X02
        X2d = X20 - X02

        Xeq_00 = templc
        Xeq_10 = 0.0
        Xeq_01 = 0.0
        Xeq_20 = cs2 * templc
        Xeq_02 = cs2 * templc
        Xeq_11 = 0.0
        Xeq_21 = 0.0
        Xeq_12 = 0.0
        Xeq_22 = cs2 * cs2 * templc
                
        Xeq_2s = Xeq_20 + Xeq_02
        Xeq_2d = Xeq_20 - Xeq_02

        X00 = X00 + wt0 * (Xeq_00 - X00)
        X10 = X10 + wt1 * (Xeq_10 - X10)
        X01 = X01 + wt2 * (Xeq_01 - X01)
        X2s = X2s + wt3 * (Xeq_2s - X2s)
        X2d = X2d + wt4 * (Xeq_2d - X2d)
        X11 = X11 + wt5 * (Xeq_11 - X11)
        X21 = X21 + wt6 * (Xeq_21 - X21)
        X12 = X12 + wt7 * (Xeq_12 - X12)
        X22 = X22 + wt8 * (Xeq_22 - X22)

        X20 = 0.5 * (X2s + X2d)
        X02 = 0.5 * (X2s - X2d)

        Xp00 = X00
        Xp10 = X10 + uilc * X00
        Xp01 = X01 + ujlc * X00
        Xp20 = X20 + uilctwo * X10 + uilc2 * X00
        Xp02 = X02 + ujlctwo * X01 + ujlc2 * X00
        Xp11 = X11 + uilc * X01 + ujlc * X10 + uilc * ujlc * X00
        Xp21 = X21 + uilctwo * X11 + ujlc * X20 + uilc2 * X01 + uilctwo * ujlc * X10 + uilc2 * ujlc * X00
        Xp12 = X12 + ujlctwo * X11 + uilc * X02 + uilctwo * ujlc * X01 + ujlc2 * X10 + uilc * ujlc2 * X00
        Xp22 = X22 + uilctwo * X12 + ujlctwo * X21 + 4.0 * uiujlc * X11 + ujlc2 * X20 + uilc2 * X02 + uilc2 * ujlctwo * X01 + uilctwo * ujlc2 * X10 + uilc2 * ujlc2 * X00

    
        t2[ina0] = Xp00 - Xp20 - Xp02 + Xp22
        t2[ina1] = 0.5 * (Xp10 + Xp20 - Xp12 - Xp22)
        t2[ina2] = 0.5 * (Xp01 + Xp02 - Xp21 - Xp22)
        t2[ina3] = 0.5 * (-Xp10 + Xp20 + Xp12 - Xp22)
        t2[ina4] = 0.5 * (-Xp01 + Xp02 + Xp21 - Xp22)
        t2[ina5] = 0.25 * (Xp11 + Xp21 + Xp12 + Xp22)
        t2[ina6] = 0.25 * (-Xp11 + Xp21 - Xp12 + Xp22)
        t2[ina7] = 0.25 * (Xp11 - Xp21 - Xp12 + Xp22)
        t2[ina8] = 0.25 * (-Xp11 - Xp21 + Xp12 + Xp22)
        end
    end    
end