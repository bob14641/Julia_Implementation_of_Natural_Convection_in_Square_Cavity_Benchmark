function vout(nnx,nny,solidstr,rho,rhoold,ui,uiold,uj,ujold,temp,tempold,n)

g2r = 0.0
g2ui = 0.0
g2th = 0.0
g2uj = 0.0
cnt = 0

for x = 1 : nnx
    for y = 1 : nny
        in = indicies_scalar(nnx,x,y)
        if solidstr[in] == 0
            g2r += (rho[in] - rhoold[in]) * (rho[in] - rhoold[in])
            g2ui += (ui[in] - uiold[in]) * (ui[in] - uiold[in])
            g2uj += (uj[in] - ujold[in]) * (uj[in] - ujold[in])
            g2th += (temp[in] - tempold[in]) * (temp[in] - tempold[in])

            cnt += 1

            rhoold[in] = rho[in]
            uiold[in] = ui[in]
            ujold[in] = uj[in]
            tempold[in] = temp[in]
        end
    end
end

g2r = sqrt(g2r / cnt)
g2ui = sqrt(g2ui / cnt)
g2uj = sqrt(g2uj / cnt)
g2th = sqrt(g2th / cnt)

global errorx = g2ui
global errory = g2uj
global errort = g2th
global errorr = g2r

if rem(n,1000) == 0
    println("nstep ", n ,"|| g2th ", g2th ,"|| g2ui ",g2ui ,"|| g2uj ",g2uj ,"|| g2r ",g2r)
end

end