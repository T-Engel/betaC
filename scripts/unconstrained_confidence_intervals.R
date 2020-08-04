alpha_k_m= function(k,n,m){
    if(k<=n-m){out= choose(n-k, m)/choose(n,m)
    }else{out=0}
    return(out)
}

S_ind=function(x,m){
    n=sum(x)
    counts<-as.data.frame(table(x))
    k=as.numeric(as.character(counts[,1]))
    f_k= as.numeric(counts[,2])

    alphas=sapply(k, alpha_k_m, n=n,m=m)

    out= vegan::specnumber(x)-sum(alphas*f_k)
    return(out)

}

var_rare= function(x,m){
    counts<-as.data.frame(table(x))
    k=as.numeric(as.character(counts[,1]))
    f_k= as.numeric(counts[,2])
    alphas=sapply(k, alpha_k_m, n=n,m=m)

    sum((1-alphas)^2*f_k)- S_ind(x,m=m)^2/calc_chao1(x)




}
c(S_ind(x,21), S_ind(x,21)-1.96*var_rare(x,21))


