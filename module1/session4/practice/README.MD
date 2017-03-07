> **Author**: Adolfo Hoyos Díaz
>
> Pontificia Universidad Javeriana, Cali



# Practice: BIOCHAM 

**Advisor**: Diana Hermith Ph.D.

## 1. Reversible enzyme kinetics



The kinetics for this example are described by the equations:



$S+E \array{ k_1 \\ \longleftrightarrow  \\ k_2 } ES$

$ES \array{ k_3 \\ \longleftrightarrow  \\ k_4 } EP $



$EP \array{ k_5 \\ \longleftrightarrow  \\ k_6 } P+E $



The designed code for the proposed dynamics is:

```c
MA (K1) for S + E => ES.
MA (K2) for ES => E + S.
MA (K3) for ES => EP.
MA (K4) for EP => ES.
MA (K5) for EP => P + E.
MA (K6) for P + E => EP.
parameter(K1, 1).
parameter(K2, 1).
parameter(K3, 1).
parameter(K4, 1).
parameter(K5, 1).
parameter(K6, 1).
present(S, 1).
present(E, 1).
```
And put inside Biocham:

![][rever]

### Simulations

#### ODE 

When the code is run with ODE simulator, using the Rosenbrok algorithm, it generates a plot like the one shown next:

![][reverrosen]

Now, selecting the Runge-Kutta algorithm for solving the ODE system: 

![][reverrunge]

Although it is not directly visible, the generated graphics through both methods almost exactly overlap. This tells us that, at this level of complexity and scale, both methods behave practically the same.

#### Stocastic



The next plot was generated selecting the Gillespie algorithm:

![][reveGillespie]



On the other hand, the Tau-Lipping algorithm was selected to get the next figure:

![][reveTau]



While the results of ODE and stocastic simulations are apparently different, both Gillespie and Tau-Lipping curves oscillate around the obtained ODE plots, and present similar steady states for the different observable variables.

#### Boolean

The boolean simulation was ran during 30 transitions, and its plot is shown on the next figure:

![][reveBoolean]

## 2. Coupled irreversible enzyme kinetics



The reaction equations are:

$S+E \array{ k_1 \\ \longleftrightarrow  \\ k_2 } ES$

$ES \array{ k_3 \\ \longrightarrow  } P+E $

$P+E \array{ k_5 \\ \longleftrightarrow  \\ k_4 } EP $

$EP \array{ k_6 \\ \longrightarrow  } S+E $

And the respective Biocham code is:

```
MA (K1) for S + E => ES.
MA (K2) for ES => E + S.
MA (K3) for ES => P+E.
MA (K4) for EP => P+E.
MA (K5) for P + E => EP.
MA (K6) for EP => S + E.
parameter(K1, 1).
parameter(K2, 1).
parameter(K3, 1).
parameter(K4, 1).
parameter(K5, 1).
parameter(K6, 1).
present(S, 1).
present(E, 1).
```



![][irrever]

### Simulations



#### ODE

![][irreODE]

#### Stocastic



![][irreGillespie]







![][irreTau]

In a similar way than the *reversible enzyme kinetics* example, all curves in both stocastic simulations, for the *irreversible enzyme kinetics* example, oscillate aroung the ODE plot.



#### Boolean

The boolean simulation is shown next:

![][irreBoolean]



## Comparisson

It can be seen, with the values of all parameters set to "1", that the shape of the dynamics of all components are very similar between both the reversible and the irreversible dynamics. The only apparent differences in this case is the stabilization time, which is less for the second case. This can be explained due to the irreversible dynamics proposed for the reaction in which we can see two irreversibles stages governed by the parameters $k_3$ and $k_6$ :

$\array{ S+E &  \array{ k_1 \\ \longleftrightarrow  \\ k_2 } & ES \\   \array{ \uparrow  & k_6} & & \array{ k_3  & \downarrow  } \\ EP & \array{ k_4 \\ \longleftrightarrow  \\ k_5 } & P+E }  $                        



In contrast, the reversible dynamics is described only by two-way reactions:



$  \array{ S+E   &   \array{ k_1 \\ \longleftrightarrow  \\ k_2 }  & ES  \\   &  & \array{ k_3 & \updownarrow & k_4} \\  EP &  \array{ k_6 \\ \longleftrightarrow  \\ k_5 }  & EP  }$



More analytical and numerical tests have to be done to precisely address the causes of this differences, e.g. with different parameter values.


[rever]: images/biocham_reversible1.PNG
[reverrosen]: images/reversibleRosenbrok1.png
[reverrunge]: images/reversibleRunge1.png
[reveGillespie]: images/reveGillespie.png
[reveTau]: images/reveTau.png
[reveBoolean]: images/reveBoolean.png



[irrever]: images/biocham_irreversible1.PNG
[irreODE]: images/irreODE.png
[irreGillespie]: images/irreGillespie.png
[irreTau]: images/irreTau.png
[irreBoolean]: images/irreBoolean.png