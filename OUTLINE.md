List of coding tasks to implement between now, 4/22/2020, to Monday 4/27/2020.

Be able to model algebraically the rate in which a disease spreads.
Be able to then add more detail to this model (still algebraically)
  -speed in which individuals travel
  -quarantine
  -population density
We are in hopes to finish the general model by Monday.



I am currently trying to figure out how to integrate the functions defined by the SIR model in to the rootsolver.py function. I'm thinking of placing a bunch of if/else statements in the functions themselves to have different values based on the time t. The main issue I'm running in to is that it appears the function doesn't change signs, making rootfind not ideal. I think that I well end up having to rip code from last weeks session. 

After checking, the rootfind.py function doesn't really do what I need it to. Switching over to the code that I wrote for week 9. 