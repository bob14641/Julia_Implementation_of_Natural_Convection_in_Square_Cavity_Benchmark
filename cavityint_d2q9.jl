function cavityint_d2q9(nnx,nny,rhoo,rho,rhoold,ui,uiold,uj,ujold,pres,temp,tempold,tmid,fri,frj,solidnode,solidstr,cs2)

for x = 1 : nnx
    for y = 1 : nny
        in = indicies_scalar(nnx,x,y)
        ui[in] = 0.0
        uj[in] = 0.0
        rho[in] = rhoo
        pres[in] = cs2 * rho[in]
        temp[in] = tmid
        tempold[in] = temp[in]
        rhoold[in] = rho[in]
        uiold[in] = ui[in]
        ujold[in] = uj[in]
        fri[in] = 0.0
        frj[in] = 0.0

        solidstr[in] = 0
        solidnode[in] = 0 
        #Bottom Wall
        if x == 1
            solidnode[in] = 1
        end
        #Left Wall
        if y == 1 
            solidnode[in] = 2
        end
        #Top Wall
        if x == nnx
            solidnode[in] = 3
        end
        #Right Wall
        if y == nny
            solidnode[in] = 4
        end
        #Bottom-Left Corner
        if x == 1 && y == 1
            solidnode[in] = 5
        end
        #Bottom-Right Corner
        if x == 1 && y == nny
            solidnode[in] = 6
        end
        #Top-Left Corner
        if x == nnx && y == 1
            solidnode[in] = 7
        end
        #Top-Right Corner
        if x == nnx && y == nny
            solidnode[in] = 8
        end
    end
end

end

