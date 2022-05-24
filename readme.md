[cover](./cover.jpg)

# MACS Adaptive Regulation Benchmark (MACS ARB)

https://github.com/macs-lab/macs-arb-adaptive-regulation-benchmark

The MACS ARB is a top-tier algorithm in the 2012-2013 [*International Benchmark on Adaptive Regulation*](http://dx.doi.org/10.1016/j.ejcon.2013.05.007) by Ioan D. Landau. The benchmark provided state-of-the-art evaluation and dissemination of adaptation methods for active vibration control and noise control. 

Algorithmically, MACS ARB presents an adaptive control scheme for identifying and rejecting unknown and/or time-varying narrow-band vibrations. We discuss an idea of safely and adaptively inverting a (possibly non-minimum phase) plant dynamics at selected frequency regions, so that structured disturbances (especially vibrations) can be estimated and canceled from the control perspective. By taking advantage of the disturbance model in the design of special infinite-impulse-response (IIR) filters, we can reduce the adaptation to identify the minimum amount of parameters, achieve accurate parameter estimation under noisy environments, and flexibly reject the narrow-band disturbances with clear tuning intuitions. 

See details in: [Chen, X.; and Tomizuka, M. Selective Model Inversion and Adaptive Disturbance Observer for Time-varying Vibration Rejection on an Active-suspension Benchmark. European Journal of Control, 19(4): 300 - 312. 2013. ](https://www.researchgate.net/publication/259137172_Selective_model_inversion_and_adaptive_disturbance_observer_for_rejection_of_time-varying_vibrations_on_an_active_suspension)

## Running the Code

run standard_test_embedFun_submit.m and configure different options from the command prompt.

