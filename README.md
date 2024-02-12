This is a correction code to the article "The effect of step size on straight-line orientation" by Lana Khalidy, Orit Peleg, Claudia Tocco, L. Mahadevan, Marcus Byrne, and Marie Dacke, published in the Journal of The Royal Society Interface in 2019 (DOI: 10.1098/rsif.2019.0181). 

insect_nav.m sweeps over values for kappa of BRW turning angle distribution, and weight of CRW component (w). Step size and, arena size, and CRW kappa are fixed, and determined from the experiments.

Inputs:
1) BRW_array: vector of BRW kappa values
2) w_array: vector of weight values
3) CRW_array: vector of CRW kappa values
4) step_size_array: vector of step sizes
5) N_iter: number of random walks to simulate

Output:
1) R_vec_length_mean: array of mean resultant lengths

