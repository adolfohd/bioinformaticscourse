# Module 1: Final session

This session brings a mini project described in [this file].

# Answers

## 1. 


The **exact** dynamics has the next reactions

> **ES bind** reaction: 
>
> $E + S \longleftrightarrow ES$
>
detailed as:
> $E(s) + S(Y \~{} 0) \array{k_{p_1} \\   \longleftrightarrow \\ k_{m_1} } E(s!1).S(Y  \~{} 0!1)$ 
 $k_{p_1} = \frac {1} {N_{Av}}$
 $k_{m_1}=\frac{1}{10}$


 > **ESconvert**: 
 >
 > $ES \rightarrow E+S$
 >
 > detailed as:
 > $\array{ \\E(s!1).S(Y~0!1)} \array{k_2 \\ \longrightarrow} \array{ \\E(s) + S(Y \~{}P) }$
  $k_{2}= \frac{1}{100}$

With the initial concentrations

> $ E_0 =\frac{1}{100}\times N_{Av}$
$S_0 = 1\times N_{Av}$ 

Meanwhile the **approx** file describes only one rection:
> $E+S\rightarrow E+S$
>
> detailed as:
> $\array{ \\ E(s) + S(Y \~{} 0)} \array{\frac{kcat}{(Km + Su)} \\ \longrightarrow} \array{ \\E(s) + S(Y \~{} P)  }$

with the same initial concentrations.

We can see that both reactions start with the same concentrations of enzyme and **dephosporylated** substrate and, although the *exact reaction* describes an intermediate step in which the initial seed converts into the **ES** substance, both reactions end with enzyme and phosphorylated substrate. We are yet to check, through simulation, whether the result is the same or not... most probably not:

The reaction rate of the **approx** file

$MM_{rate} =\frac{kcat}{(Km + Su)}$

depends on the concentration of dephosporylated compound SU; this will speed up the reaction since the rate will grow as the concencentration of SU decreases.


| Feature |Units|Exact | Aprox |
| --- | --- | --- | --- |
| Volume | ${um}^3$ | 1 | 1


### Simulation results

#### Exact

![exact][]

We can see that the transition state is very ephimeral, so the concentration of ES compound remains very low during the whole process.
 
#### Approx

![approx][]

[this file]: project.pdf
[exact]: 1/exact.png
[approx]: 1/approx.png