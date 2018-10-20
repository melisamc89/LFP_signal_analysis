function mi=ModularityIndex(phase,amplitude)

    h=-sum((amplitude+eps).*log(amplitude+eps));
    
    n=length(phase);
    mi=(log(n)-h)/log(n);

end