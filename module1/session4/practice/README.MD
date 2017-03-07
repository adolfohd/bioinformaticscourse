# Practice: BIOCHAM 



## 1. Reversible enzyme kinetics

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

When the code is run with ODE using the Rosenbrok algorithm, it generates a plot like the one shown next:

![][reverrosen]

Now, selecting the Runge-Kutta algorithm for solving the ODE system: 

![][reverrunge]

Although it is not directly visible, the generated graphics through both methods almost exactly overlap. This tells us that, at this level of complexity and scale, both methods behave practically the same.

## 2. Coupled irreversible enzyme kinetics



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


[rever]: images/biocham_reversible1.PNG
[reverrosen]: images/reversibleRosenbrok1.png
[reverrunge]: images/reversibleRunge1.png
[irrever]: images/biocham_irreversible1.PNG