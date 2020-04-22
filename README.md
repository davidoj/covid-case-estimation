Estimate true numbers of cases from Johns Hopkins data.

The "ascertainment" function does most of the calculation.

## Align confirmed cases

The first step is to calculate new infections timeshifted back by 7 days. The figure 7 comes from:
 - 5.5 days infection-to-symptoms (https://docs.google.com/spreadsheets/d/1yzVSp71jiCsoD_L6sXchfg8L9OF0tGHWEnCdfvF73Ac/edit#gid=1242721729)
 - 4 days symptoms-to-confirmation (also https://docs.google.com/spreadsheets/d/1yzVSp71jiCsoD_L6sXchfg8L9OF0tGHWEnCdfvF73Ac/edit#gid=1242721729)
 - Subtract two days as exponential growth means that fast confirmations will be overrepresented
 
We could in principle calculate this properly as we do with time to death, but for the purposes of calculating the true number of cases it is less crucial than calculating infection to death properly.

There may also be a lag in the reporting of deaths, and this would have a significant impact on results.

## Comparison of Expected Deaths

The number of deaths on day $i$ should be
$$ d_i = f \sum_j p_{i-j} t_j $$
where $f$ is the infection fatality rate, $p_k$ is the probability of dieing $k$ days after being infected and $t_j$ is the true number of cases on day $j$.

Suppose that there is a fixed ascertainment rate $a$ such that confirmed cases $c_j=a t_{j-7}$ for all $j$. Then 

$$ a = \frac{f\sum_j p_{i-j} c_{j+7}}{d_{i}} $$

**Note** This is actually *lagged* ascertainment rate - it asks "what fraction of new cases today are detected 7 days later". I don't directly estimate an instantaneous ascertainment rate as this is very sensitive to the short term trajectory of infections, which, as I note above, I can't yet forecast. What we can in fact calculate is:

$$ t_{i-7} = \frac{c_i}{a}  $$


# Filter out small numbers of deaths

The method doesn't work well with small numbers of deaths per day. There are probably better ways to do this like variance weighted averaging, but I just remove small numbers of deaths.


## Recovered, Infectious and Exposed Classes

This is for use if you are trying to determine initial conditions for an SEIR model.

The infectious class is given by 

$$ I_i = \sum_{j} q_{i-j} t_j - \sum_{j} r_{i-j} t_j $$

and the exposed class by 

$$ E_i = \sum_j (1-q_{i-j}) t_j $$

finally, the recovered class is given by

$$ R_i = \sum_j r_{i-j} t_j $$

Where $q_n$ is the distribution of the time from being exposed to being infectious, $E_i$ and $I_i$ the size of the exposed and infected classes on day $i$ and $t_j$ the true number of new cases on day $j$ as above

To make the infectious class match standard SEIR treatment, we take $q_{i-j}=A e^{-(i-j)\epsilon}$ where $A$ is some normalising constant. That is, the transition probability is an exponential. Currently we use $\epsilon=0.2$.

To get a recovery rate that matches our use of GLEAM, we could take $r_i = B e^{-i\mu}$ where $\mu=0.5$ is the rate of recovery from GLEAM and $B$ is a normalising constant. This is a bit unrealistic, but alternative approaches are tricky to implement.

Another possibility is to take the "serial interval" distribution from https://www.imperial.ac.uk/media/imperial-college/medicine/mrc-gida/2020-03-30-COVID19-Report-13.pdf (note that this seems to be the distribution of probability of infection, not probability of first infection, so it's not acutally the serial interval). However, we need to turn this nonexponential distribution into an exponential infectious class, and figuring out how to approximate this appropriately seems hard.