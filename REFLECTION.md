# Reflection

## The Model
The model state consists of the current `x`,`y` position of the car, its angle `psi`, its velocity `v` as well as the cross track error `cte` and the heading error `epsi` (MPC.cpp: 137-142).

The update equations are taken from the lecture and are implemented at lines 115-120 (MPC.cpp).

The optimization finds the best values for the actuators `steering_angle` and `throttle`, whose values are used in the update equation.

### Error weights
Most time for this project was spent on finding proper weights for the different error terms.
The most important of these are the `cte` and `epsi` weights, which are both set to `20`.
The speed deviation weight is set at `1`, while the remaining weights for the minimization of steering and throttle as well as for minimizing the gaps between consecutive values of these actuators are all set to `0.1`.
These values are a result of a trial-and-error approach, where the error and actuator values were plotted and weights were sought that minimized both the error and oscillatory behaviour (main.cpp: 194).

This was done by making the weights configurable from the command line (main.cpp: 79:97).

## Timestep length and elapsed duration
The parameters `N` and `dt` were chosen as 10 and 0.1, respectively.
This choice was influenced by the latency of 100ms.
Predicting for one second into the future then followed naturally.
These values seem to be good ones, as varying them did not improve the controller markedly.

## Polynomial Fitting and MPC Preprocessing
The waypoints in `ptsx` and `ptsy` come in the global coordinate system of the simulator.
They are therefore transformed into the coordinate system of the car (main.cpp: 150-162) and visualized as the yellow reference line.

A polynomial of grade 3 is then fitted to the waypoints and used in the MPC solution.

The resulting steering value is multiplied by -1 and scaled to +/- 25 degrees (in radians).

## Model Predictive Control with Latency
The latency is modelled by using a future location of the car as the input to the MPC (main.cpp: 135-148).

The position, heading and velocity of the car are updated by predicting them with an offset of half the latency (main.cpp: 143).
This produced much stabler results than using the full latency during the prediction.


# Simulation
The car successfully completes laps at around 75Mph, as can be seen in the video `result.ogv`.
