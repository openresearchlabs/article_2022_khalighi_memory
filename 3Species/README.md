# Reproduce the dynamics of three-species microbial communities with arbitrary conditions

## 3specis

 ### README.md              # project overview

 ### Main_ThreeSpecies.m      # Run this code

 *  **Method1.m**            # Numerical methods for the model under given perturbations
 * **FDE_PI12_PC.m**        # solver for fractional differential equation  
 *  **OrnsteinUhlenbeck.m**  # Ornstein Uhlenbeck stochastic process $dX_t = \theta (\mu - X_t)dt + \sigma dW_t$
 *  **b3OUP.mat**            # Growth rate generated under stochastic pertubation used in the paper, 700x3
