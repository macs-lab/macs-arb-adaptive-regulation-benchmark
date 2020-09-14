%% closed loop sys test
if SW_BASELINE_CONTROL_SYS
    Gc      = tf(R,S,Ts);
    Gp      = tf(B,A,Ts);
    HOL     = Gc*Gp;
    HCL     = minreal(HOL/(1+HOL));
    
    figure, bodeplot(HOL,bode_opt)
    grid on;zoom on;
    title 'Baseline open loop frequency response',xlim([0,400])
    
    figure, bodeplot(HCL,bode_opt)
    grid on;zoom on;
    title 'Baseline closed loop frequency response',xlim([0,400])
    
    figure, bodeplot(tf(C,D,Ts),bode_opt)
    grid on;zoom on;
    title('Frequency response of the primary path'),xlim([0 400]);
    
    figure, step(HCL),  title('Baseline closed loop step response')
    figure, step(HOL),  title('Baseline open loop step response')
    figure, pzplot(HCL),title('Baseline closed loop pole zero plot')
    figure, pzplot(HOL),title('Baseline open loop pole zero plot')
end