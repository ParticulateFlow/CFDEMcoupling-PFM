{
    m=gSum(rhoeps*1.0*rhoeps.mesh().V());
    if(counter==0) m0=m;
    counter++;
    Info << "\ncurrent gas mass = " << m << "\n" << endl;
    Info << "\ncurrent added gas mass = " << m-m0 << "\n" << endl;
}
