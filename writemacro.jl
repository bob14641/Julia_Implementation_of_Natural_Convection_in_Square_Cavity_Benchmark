function writemacro(nnx,nny,ui,uj,rho,temp,printui,printuj,printrho,printtemp,solidnode,solidnodea,n)

    for x = 1 : nnx
        for y = 1 : nny 
            in = indicies_scalar(nnx,x,y)
            printui[x,y] = ui[in] 
            printuj[x,y] = uj[in]
            printrho[x,y] = rho[in]
            printtemp[x,y] = temp[in]
            solidnodea[x,y] = solidnode[in]
        end
    end

    if n == 1 
        debuga=open("julia_m_solidnode_debug.txt","w")
            writedlm(debuga, solidnodea, ",")
        close(debuga)
        xvelo = open("julia_m_xvelo_$n.txt", "w") 
            writedlm(xvelo, printui , ",")
        close(xvelo) 
        yvelo = open("julia_m_yvelo_$n.txt", "w") 
            writedlm(yvelo, printuj , ",")
        close(yvelo)
        rhoa=open("julia_m_rho_$n.txt","w")
            writedlm(rhoa, printrho, ",")
        close(rhoa)
        tempa=open("julia_m_temp_$n.txt","w")
            writedlm(tempa, printtemp, ",")
        close(tempa)
    else 
        errn = open("julia_m_$n.txt", "w") 
            write(errn, "n = $n, g2th = $errort, g2ui = $errorx, g2uj = $errory, g2r = $errorr \n")
        close(errn)
        xvelo = open("julia_m_xvelo_$n.txt", "w") 
            writedlm(xvelo, printui , ",")
        close(xvelo) 
        yvelo = open("julia_m_yvelo_$n.txt", "w") 
            writedlm(yvelo, printuj , ",")
        close(yvelo)
        rhoa=open("julia_m_rho_$n.txt","w")
            writedlm(rhoa, printrho, ",")
        close(rhoa)
        tempa=open("julia_m_temp_$n.txt","w")
            writedlm(tempa, printtemp, ",")
        close(tempa)
    end
end