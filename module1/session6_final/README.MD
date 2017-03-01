# Module 1: Final session

This session brings a mini project described in [this file].

# Answers

## 1. Exact vs Approx.


The **exact** dynamics has the next reactions

> **ES bind** reaction: 
>
> $E + S \longleftrightarrow ES$
>
detailed as:
> $E(s) + S(Y \sim 0) \array{k_{p_1} \\   \longleftrightarrow \\ k_{m_1} } E(s!1).S(Y  \sim 0!1)$ 
with
>
 $k_{p_1} = \frac {1} {N_{Av}}$
 $k_{m_1}=\frac{1}{10}$

>
>
 > **ESconvert**: 
 >
 > $ES \rightarrow E+S$
 >
 > detailed as:
 > $\array{ \\E(s!1).S(Y \sim 0!1)} \array{k_2 \\ \longrightarrow} \array{ \\E(s) + S(Y \sim P) }$
 >
 with
 >
  $k_{2}= \frac{1}{100}$

With the initial concentrations

> $ E_0 =\frac{1}{100}\times N_{Av}$
$S_0 = 1\times N_{Av}$ 

Meanwhile the **approx** file describes only one rection:
> $E+S\rightarrow E+S$
>
> detailed as:
> $\array{ \\ E(s) + S(Y \sim 0)} \array{\frac{kcat}{(Km + Su)} \\ \longrightarrow} \array{ \\E(s) + S(Y \sim P)  }$
>
>with
>
$k_f=\frac{1}{N_{Av}}$
$k_r=\frac{1}{10}$
$k_{cat}=\frac{1}{100}$
$k_m=\frac{k_{cat}+k_r}{k_f}$


with the same initial concentrations.

We can see that both reactions start with the same concentrations of enzyme and **dephosporylated** substrate and, although the *exact reaction* describes an intermediate step in which the initial seed converts into the **ES** substance, both reactions end with enzyme and phosphorylated substrate. We are yet to check, through simulation, whether the result is the same or not... most probably not:

The reaction rate of the **approx** file

$MM_{rate} =\frac{kcat}{(k_m + Su)}$

depends on the concentration of dephosporylated compound SU; this will speed up the reaction since the rate will grow as the concencentration of SU decreases.


| Feature |Units|Exact | Aprox |
| --- | --- | --- | --- |
| Volume | ${um}^3$ | 1 | 1


### Simulation results

#### Exact

![exact][]

Stability of the process is reached around $1.4\times 10 ^{4} s$, and also we can see that the transition state is very ephimeral, so the concentration of ES compound remains very low during the whole process, as shown on the next plot:

![exactES][]

ES compund starts from $\sim5.5$ and goes to $0$.


#### Approx

despite it may not ve obvious, the next graphic, obtained simulating the **approx** file, (almost) exactly overlaps with the previous graphic.

![approx][]


We apreciate that the ES compound is being monitored as an observable variable, even when this substance is not present in the described reaction. This leads to the next graphic:

![approxES][]

that shows this variable always as zero, as should be expected.

The results for the $MM_{rate}$ function is seen in the next figure

![approxMMrate][] 
 

and, although the wave form relates (almost) exactly to the ES evolution of the **exact** reaction, the values are two orders of magnitud below, going from $\sim 5.5 \times 10 ^{-2}$ to $0$.


## 2. A little modification

For this part, we added the function:

> $ES_{eff} = E_0 \times \frac{Su}{k_m + S_u}$

let's remember we have

> $k_f=\frac{1}{N_{Av}}$
$k_r=\frac{1}{10}$
$k_{cat}=\frac{1}{100}$
$k_m=\frac{k_{cat}+k_r}{k_f}$


so

> $ES_{eff} = E_0 \times \frac{Su}{\frac{k_{cat}+k_r}{k_f} + Su}$
The graph obtained for this function from the simulation results is shown next, and it overlaps with the ES complex of the **exact** reaction.


![ESeff][]

![ESMM][]

[this file]: project.pdf
[exact]: 1/exact.png
[exactES]: 1/exactES.png
[approx]: 1/approx.png
[approxES]: 1/approxES.png
[approxMMrate]: 1/approxMMrate.png

[ESMM]: 2/approxESMM.png
[ESeff]: 2/approxESeff.png