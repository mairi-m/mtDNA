# Agent-based models of mtDNA population dynamics
This project is about modelling clonal expansion. It is for the [gold CREST award](https://www.crestawards.org/crest-gold) with [Conor Lawless](https://www.ncl.ac.uk/medical-sciences/people/profile/conorlawless.html) as a mentor.

## Mitochondria
[Mitochondria](https://en.wikipedia.org/wiki/Mitochondrion) are often referred to as 'the powerhouse of the cell'. Mitochondria are a type of sub-unit ([organelle](https://en.wikipedia.org/wiki/Organelle)) within a [cell](https://en.wikipedia.org/wiki/Cell_(biology)) which create energy for a cell using oxygen, fats, sugars and proteins. Mitochondria are also involved in cell death ([apoptosis](https://en.wikipedia.org/wiki/Apoptosis)) which is important because cells which divide too many times result in exessive numbers of cells, causing tumors.

![Animal cell](https://upload.wikimedia.org/wikipedia/commons/4/48/Animal_cell_structure_en.svg)

> Cells have different organelles within them, mitochondria are one of them.

![Mitochondria](images/Mitochondrion_structure.png)

> Diagram of a mitochondrion.

![mitochondrion](images/mitochondrion.webp)

> Image of a mitochondrion

Malfunctioning mitochondria can have [negative health effects](https://www.umdf.org/what-is-mitochondrial-disease-2/0) such as fatigue, muscle weakness and neurodegenerative diseases since they are responsible for making the majority of energy in the body.

## Mitochondrial DNA
Cells which have mitochondria ([eukaryotic](https://biologydictionary.net/eukaryotic-cell/)) not only have DNA in the nucleus but have mitochondrial DNA (mtDNA) in each mitochondrion. mtDNA is significantly smaller than DNA in the nucleus, however, because cells almost always have many mitochondria and a mitochondrion almost always has multiple mtDNA, there is many more copies of mtDNA in a cell than nuclear DNA. mtDNA is circular and more similar to bacterial DNA than linear nuclear DNA.

Mutated mtDNA can cause disease if they cannot carry out their function properly.

## Clonal Expansion
[Clonal expansion](https://royalsocietypublishing.org/doi/10.1098/rsob.200061) is when the concentration mtDNA mutations within cells increase. It can happen as humans age or because of mitochondrial disease. Larger amounts of mutated mtDNA can also cause mitochondrial diseases.

## Mitochondrial Disease
Mutated mtDNA molecules can cause disease because they code for many mitochondrial proteins.  If mutant mtDNA species reaches sufficiently high levels in the cell then these proteins will not be produced in sufficient quantities to support cell function. The result is mitochondrial dysfunction. Some primary mitochondrial diseases are [leigh syndrome](https://en.wikipedia.org/wiki/Leigh_syndrome), [mitochondrial myopathy](https://en.wikipedia.org/wiki/Mitochondrial_myopathy) and [Leber's hereditary optic neuropathy](https://en.wikipedia.org/wiki/Leber%27s_hereditary_optic_neuropathy). Mitochondrial dysfunction has also been connected to illness such as alzheimers, cancer and muscular dystrophy. Mitochondrial disease can have severe and untreatable consequences such as vision and hearing loss, difficulty walking and breathing, epilepsy, neurodegeneration and developmental defects.

## Modelling mtDNA
It is important to model mtDNA to gain a better understanding of it so treatment for many different mitochondrial diseases can be developed. In this project we will create a model of skeletal muscle fiber with wild type and mutant mtDNA diffusing randomly. The model will show mtDNA being synthesised and degraded. It will also show perinuclei which are fixed in space and how the mtDNA interacts with them.

![skeletal muscle fiber](images/skeletal_muscle_fiber_d.jpg)
> mtDNA travels slowly longitudinally and quickly radially.

## Tools
We will create [an agent based model](https://en.wikipedia.org/wiki/Agent-based_model) using [python](https://www.python.org/) and the [MESA library](https://mesa.readthedocs.io/en/stable/). Python is a popular easy to read programming language which I already have experiance in. It can also be used with the MESA library for agent based modelling which is why we are using it.

## Timescale

In the beginning I split the project into 5 sections:
- Meetings with Conor Lawless
- Researching the biology behind the project
- Learning how to use the mesa library
- Creating the model
- Documentation

I researched the biology first to gain a better understanding of the project and then documented what I learned in github. I then began learning how to use the mesa library. I had meetings with Conor and added to the documentation throughout the project. 

I used google sheets ([sheet linked here](https://docs.google.com/spreadsheets/d/1H-vyV8cI5vYHMRJBZxGlX5kHtQQPEzIrgFyDWQjZb1U/edit?usp=sharing)) to organise my time and ensure I spent long enough on the project.



## References
http://mito.ncl.ac.uk/clonexp/clonal_expansion/




