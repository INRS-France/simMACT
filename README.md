# simMACT
Une bibliothèque logicielle de calcul des couples d'actionnement articulaires humains maximaux pour l'aide à la conception de poste de travail ergonomiques.
*A software library to compute human joint actuation torques for ergonomics purpose (english version below).*

`simMACT` (simulation of Maximal ACtuation Torques) est un projet de bibliothèque logicielle issue des travaux de recherche menés par l'INRS et Inria-AUCTUS. Ces travaux concernent l'estimation des performances physiques humaines pour la conception d'équipements de travail. L'objectif de `simMACT` est d'améliorer les fonctionnalités de simulation des efforts articulaires d'un logiciel de mannequin numérique (Digital Human Model, DHM) en y associant, comme une boite noire, un modèle biomécanique avancé.

## Architecture du projet
`simMACT` repose sur deux éléments de base:
- le formalisme mathématique des zonotopes, une classe particulière de polytopes, efficace pour modéliser l'actionnement d'un système rigide poly-articulé;
- un moteur de simulation musculo-squelettique.

À titre de démonstrateur, `simMACT` utilise le moteur de simulation musculo-squelettique **[OpenSim](https://opensim.stanford.edu/)** à partir de sa version 4.0. Le projet `simMACT` est implémenté en langage Python 3.8 pour des raisons de compatibilité avec cette version de **OpenSim**. 

## Diffusion et licence
`simMACT` sera prochainement disponible depuis ce site, selon aux termes de la licence BSD-3 clauses.

# simMACT
`simMACT` (simulation of Maximal ACtuation Torques) is a software library resulting from research carried out by INRS and Inria-AUCTUS. This work concerns the estimation of human physical performance for the design of work equipment. The aim of `simMACT` is to improve some simulation features of a digital mannequin software (Digital Human Model, DHM), namely the computation of joint actuation torques. It relies on an advanced biomechanical model.

## Architecture du projet
`simMACT` is based on two fundamental parts:
- the mathematical formalism of zonotopes, a particular class of polytopes, effective for modeling the actuation of a rigid polyarticulated system;
- a musculoskeletal simulation engine.

As a demonstrator, `simMACT` uses the musculoskeletal simulation engine **[OpenSim](https://opensim.stanford.edu/)** from version 4.0. The `simMACT` project is implemented in Python 3.8 for compatibility with this version of **OpenSim**.

##Distribution and license
`simMACT` will soon be available from this site under the terms of the BSD-3 clauses license.


