#!/bin/bash


#concatenate all SDs called from 1000GP reference guided assemblies
cat *.bed > SDs_from_all_samples.bed

#deduplicate SDs
sort SDs_from_all_samples.bed | uniq > HGSVC3_SDs_from_65_individuals.bed