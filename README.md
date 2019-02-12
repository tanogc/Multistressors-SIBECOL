# Modelling multi-stressor responses across ecosystems in R

Este taller es una introducción al modelado de los efectos de estresores múltiples en los ecosistemas, que tuvo lugar durante el congreso SIBECOL2019. Puedes consultar el guion en [pdf](https://github.com/tanogc/Multistressors-SIBECOL/blob/master/Guion%20del%20curso.pdf).

En este taller mostramos cómo analizar respuestas ecológicas frente a dos o más variables explicativas que puedan tener efectos combinados. Para ello usamos unan serie de técnicas estadísticas a través de un protocolo que consta de tres partes:
1. Preparación de los datos (transformación, estandarización, colinealidad)
2. Exploración de los datos (correlaciones, Random Forest, Boosted Regression Trees)
3. Modelos finales (GLM). 

También incluye una sección de interpretación de resultados y efectos combinados de los estresores.

## Archivos

El archivo 'cookbook.R' contiene el código para ejecutar los análisis estadísticos en R. Se proporciona la función *sim.multi.str()* que permite usar datos simulados, que se encuentra en el archivo 'simul_functions.R'.

```
# Simulated dataset

n <- 100 # number of sites
ac <- 3 # accuracy SD units of error

# Single stressor hierarchy s1 > s2 > s3 = s4
# Interaction hierarchy 1:3 > 1:2 > 2:4
# s1, s2, s3, s4, s1:s2, s1:s3, s2:s4
ses<-c(5, 3, 2, 2, 2, 3, 1)

# Simulating data (additive, antagonistic, synergistic, opposing and mixed responses)
sim.multi.str(n, ses, ac, mod.type="additive", plot.int=T)->sim.set.add
sim.multi.str(n, ses, ac, mod.type="antagonistic", plot.int=T)->sim.set.ant
sim.multi.str(n, ses, ac, mod.type="synergistic", plot.int=T)->sim.set.syn
sim.multi.str(n, ses, ac, mod.type="opposing", plot.int=T)->sim.set.opo
sim.multi.str(n, ses, ac, mod.type="mixed", plot.int=F)->sim.set

sim.set$sim.dat->dat
```

## Artículo original

Este protocolo está basado en [Feld et al (2016)](https://www.sciencedirect.com/science/article/pii/S0048969716314310?via%3Dihub) 

Para citar el artículo original:

```
Feld, C.K., Segurado, P., Gutiérrez-Cánovas, C., 2016. Analysing the impact of multiple stressors 
in aquatic biomonitoring data: A 'cookbook' with applications in R. Sci. Total Environ. 573. 
```

## Dependencias
Para poder seguir el taller es necesario tener instalados las librerías 'usdm', 'randomForestSRC', 'ggRandomForests', 'gbm', 'dismo', 'MuMIn', 'lattice' y 'variancePartition'. Las podemos instalar usando los siguientes comandos:

```
install.packages("usdm")
install.packages("randomForestSRC")
install.packages("ggRandomForests")
install.packages("gbm")
install.packages("dismo")
install.packages("MuMIn")
install.packages("lattice")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("variancePartition", version = "3.8")
```

## Atribución y contacto
Si utilizas o modificas el material contenido en este repositorio, cítalo adecuadamente como:

```
Gutiérrez-Cánovas, C., Capdevilla, P. (2019) Modelling multi-stressor responses across ecosystems in R
GitHub repository, https://github.com/tanogc/Multistressors-SIBECOL.git
```
Para cualquier uso que trascienda dicha licencia o cualquier aclaración, contactar con Cayetano Gutiérrez Cánovas (cayeguti@um.es) o Pol Capdevilla Lanzaco (pol.capdevila@zoo.ox.ac.uk).




