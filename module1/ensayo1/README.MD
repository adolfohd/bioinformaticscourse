

# The Complexity of Cell-Biological Systems: A review

This essay is a review of the paper [Gabbay, Dov, et al. "The Complexity of Cell-Biological Systems."]. I'm going through the proposed paper and will highlight the, in my opinion, most relevant concepts, and the ones I find new to my experience.

## Inherent complexity

Cells are the basic building blocks of the very complex and always dynamic living systems. Cells are composed by a very large number of components, are very difficult to observe, and show non-linear dynamics.

The dynamic nature of cell-biological systems is described by the cell cycle, which is the series of steps to perform cell replication; this allows living organisms to develop and to be maintained.

Living organisms as well as cells, fabricate themselves and possess certain independency with the exterior.

The functioning of a molecule, cell or system by itself does not define its function. The function of a system depends on the numbers and kinds of neighbours it interacts with.

These interactions can be modeled through the tools that systems theory provides (switch, amplifier, oscillator, filter). Unfortunatelly, these models must be restricted to a technologically feasible context or level (size, time, population, etc.), due to the high complexity of cell-biological systems. 



## Modelling

Hence, simplified hypothesis must be formulated. Dynamical systems theory provide the tools to make these formulations and to test them, through models that, in principle, should closely describe the behaviour of the system's behaviour the are intended to model; additionally, the models must be sufficiently simple to be able to be tracked and computed, given the resources available for this tasks. This inevitably brings innacuracies to the model, so it never fully describes the system's behaviour. 

The chosen level of study to generate a model brings assumptions and choices, abstractions that bring philosophical issues:

> **Ontology**: Describes basic assumptions about things. 

>  **Epistemology**:  study of the tools available for obtaining information of things. 

The a priori knowledge of a system affects the ability to learn from it.

 The assumption that biological entities interact through feedback loops can be discussed through the tools developed by systems theory.

A system can exist in a set of states, that describe its dynamics, and can be fully measured or partially estimated from measurements. 

Some of the variables of the model of the system can have fixed values given the environment of modelling; these are parameters. The *identification* of a system is the proposition of a model's architecture and the finding of its parameters based on a priori knowledge and measurements.

It is common to propose models with architectures of differential equations, although this has been criticized, due to the different levels of systems that cannot be fully modeled at once.

## Domains of study

This brings the necessity of some conceptual closure: collapsing the exterior into some characterizations or black-boxing the interior; both are always needed at some extent. 


> **Metabolism**: Catabolism(energy storage) + anabolism (energy investment)

> **Cell signalling**: Transduction of signals, molecule transformations.

> **Gene expression**:  Transformation of gene information into gene produts.

Biology is a case of organized complexity.

> **Model organism**:  a species that has been widely studied.

## Mathematical modeling cannot be automated

Interdisiplinarity is needed and special care has to be put on the visualization of data, since the high-dimmension reduction into 2 or 3 dimmensions can hide important information.

Also, the available technological tools must be taken into account from the modelling, computations and presentation of data.

## Interactions are Intertwined

Generalizations become difficult as the knowledge of a system grows; nevertheless, we must be always aware that a trully holistic approach is never feasible. 

Additionally, having knowledge of every piece of a system does not bring certainty of the system as a whole, and moreover, the functioning of a system does not reveal the dynamics of its parts.



## Mathematical modelling: balancing compromises

Even when picking a well-defined level of study, one must keep this study to a relativelly low number of components, since the analytical and methodological tools available are limited.

Two approaches are used in computer modelling:

 - Full simulation of all parts, derived from first principles.
  - Reduction of complexity through abstraction.


The latter is the addecuate approach for the study of biological systems, since first principles cannot be obtained.


## Domains of autonomy

| Size scale | Time scale  |
|---------| ---- |
| Organism | Years |
|  Organs |  Weeks and hours |
| Cells |  Hours and minutes
| Intracellular reactions | Seconds 

Stability studied from a linear systems point of view can be inapropriate, so a non-linear perpective is better suited.

## Predictable stochasticity

Randomness can come from the inability to precisely measure a signal, of due to process that is itself random.
  
Small fluctuation around a trend can be irrelevant.

# Conclussion

The ability to precisely model a system diminishes as the system becomes more complex. Abstraction becomes necessary to achieve simplication, even when abstraction itself an be a difficult task. This brings unpredictability, meaning that no model fully describes a system.
