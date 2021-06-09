# Area-normalized scaling of ReefBudget calcification and bioerosion rates for use with CoralNet
The following code generates area-normalized calcification and bioerosion rates from ReefBudget methodologies for use with CoralNet image identification labels to support the integration of estimated carbonate production rates with CoralNet.The default rates include bioerosion to be consistent with ReefBudget methodologies, but include_bioerosion can be set to FALSE at the beginning of the script to calculate gross carbonate production rates instead.

## Estimating CoralNet calcification rates for the Indo-Pacific


Taxa-specific area-normalized calcification rates (G=kg CaCO<sub>3</sub> m<sup>-2</sup> yr<sup>-1</sup>) were calculated iteratively as G=50th percentile, G<sub>lower</sub>=25th percentile, and G<sub>upper</sub>=75th percentile of a Monte-Carlo simulation (n=10,000) using randomly selected values within the range of uncertainties for each of the taxa-specific equation terms in the below equation:

  ### <img src="https://render.githubusercontent.com/render/math?math=G=(n*cf*((c+b)*s*r%2Bi))/10">

n = number of colonies per linear meter (±95%)
	
	Source: 100/s
	
s = median colony diameter (cm) (±95%)
	
	Source: NOAA Pacific Islands coral demography data (see multiple references below)
	
cf = conversion factor accounting for open space in branching morphologies (± uncertainty) 
	
	Source: proportion of 3D space occupied by branching corals (Doszpot et al. 2019) 
	
c = calcification rate coefficient (± uncertainties) 
	
	Source: ReefBudget Indo-Pacific v1.2 (Perry et al. 2018)
	
b = bioerosion rate coefficient
	
	Source: (Microbioerosion Rate)/10; ReefBudget Indo-Pacific v1.2 (Perry et al. 2018)
	
i = calcification rate intercept (± uncertainties) 
	
	Source: ReefBudget Indo-Pacific v1.2 (Perry et al. 2018)
	
r = rugosity (± uncertainty)
	
	Source: mean (±SD) rugosity for morphology (González-Barrios and Álvarez-Filip 2018) morphology is most common from NOAA Pacific Islands coral demography data
  
10 = convert units to kg CaCO<sub>3</sub> m<sup>-2</sup> yr<sup>-1</sup>

Substitutions:

1)	Genus and morphology filled in equation terms for labels at coarser taxa resolution
2)	The mean r of all morphologies was used to fill in encrusting morphologies, assuming that encrustation occurs over mean reef structural complexity
3)	<i>Heliopora</i> branching c and i were determined by substituting extension and density data for <i>Heliopora</i> (Courtney et al., in rev) into ReefBudget <i>Pocillopora</i> branching with submassive conversion factor
4)	<i>Porites</i> foliose c and i were determined by substituting extension and density data from ReefBudget <i>Porites</i> into ReefBudget Hard Coral Foliose
5)	<i>Millepora</i> columnar c and i were determined by substituting extension and density data for ReefBudget <i>Millepora</i> into ReefBudget Hard Coral Columnar
6)  Macro- and micro- bioerosion rates from ReefBudget Indo-Pacific v1.2 were applied to all other hard substrates.

## Estimating CoralNet calcification rates for the Western Atlantic


Taxa-specific area-normalized calcification rates (G=kg CaCO<sub>3</sub> m<sup>-2</sup> yr<sup>-1</sup>) were calculated iteratively as G=50th percentile, G<sub>lower</sub>=25th percentile, and G<sub>upper</sub>=75th percentile of a Monte-Carlo simulation (n=10,0 00) using randomly selected values within the range of uncertainties for each of the taxa-specific equation terms in the below equation:

$$G=(n*cf*((c+b)*s+i))/10$$

n = number of colonies per linear meter (±95%)

	Source: 100/(s/r)
	
s = median colony surface length (cm) (±95%)

	Source: CARICOMP coral links (Dulcie 2010)
	
cf = conversion factor accounting for open space in branching morphologies (± uncertainty) 

	Source: proportion of 3D space occupied by branching corals (Doszpot et al. 2019) 
		 proportion of 3D space occupied by Acropora cervicornis (Million et al. 2021)
		 
c = calcification rate coefficient (± uncertainties) 

	Source: ReefBudget Caribbean v2 (Perry and Lange 2019)
	
b = bioerosion rate coefficient

	Source: (Microbioerosion Rate)/10; ReefBudget Caribbean v2 (Perry and Lange 2019)
	
i = calcification rate intercept (± uncertainties) 

	Source: ReefBudget Caribbean v2 (Perry and Lange 2019)
	
r = rugosity (± uncertainty)

	Source: mean (±SD) rugosity for morphology (González-Barrios and Álvarez-Filip 2018) 
	
10 = convert units to kg CaCO<sub>3</sub> m<sup>-2</sup> yr<sup>-1</sup>

Substitutions:

1)	Genus and morphology filled in equation terms for labels at coarser taxa resolution
2)	The mean r of all morphologies was used to fill in encrusting morphologies, assuming that encrustation occurs over mean reef structural complexity
3)  Microbioerosion rates from ReefBudget Caribbean v2 were applied to all hard substrates.
4)  <i>Cliona delitrix</i> bioerosion rates from ReefBudget Caribbean v2 were applied to <i>Cliona delitrix</i> and the mean (±SD) rate from ReefBudget Caribbean v2 of all Clionid sponges were applied to Clionid sponges.
5)  Mean (±SD) parrotfish bite volume from ReefBudget Caribbean v2 was applied to all bite scars.

### References:

Courtney TA, Guest JR, Edwards AJ, Dizon RM. Linear extension, skeletal density, and calcification rates of the blue coral <i>Heliopora coerulea</i>. Coral Reefs. In review.

Doszpot NE, McWilliam MJ, Pratchett MS, Hoey AS, Figueira WF. Plasticity in three-dimensional geometry of branching corals along a cross-shelf gradient. Diversity. 2019 Mar;11(3):44.

González-Barrios FJ, Álvarez-Filip L. A framework for measuring coral species-specific contribution to reef functioning in the Caribbean. Ecological Indicators. 2018 Dec 1;95:877-86.

Linton, Dulcie; Bermuda Institute of Ocean Sciences (BIOS) (2010). A unified, long-term, Caribbean-wide initiative to identity the factors responsible for sustaining mangrove wetland, seagrass meadow, and coral reef productivity, February 1993 - October 1998 (NCEI Accession 0000501). NOAA National Centers for Environmental Information. Dataset. https://accession.nodc.noaa.gov/0000501.

Million WC, O’Donnell S, Bartels E, Kenkel CD. (2021) Colony-Level 3D Photogrammetry Reveals That Total Linear Extension and Initial Growth Do Not Scale With Complex Morphological Growth in the Branching Coral, <i>Acropora cervicornis</i>. Frontiers in Marine Science. 8:384.

Perry CT, Lange I, Januchowski-Hartley FA (2018) ReefBudget Indo Pacific: online resource and methodology. Retrieved from http://geography.exeter.ac.uk/reefbudget/

Perry CT and Lange ID (2019) ReefBudget Caribbean v2: online resource and methodology. Retrieved from http://geography.exeter.ac.uk/reefbudget/

NOAA Pacific Islands Coral Demography Data:

Coral Reef Ecosystem Program; Pacific Islands Fisheries Science Center (2017). National Coral Reef Monitoring Program: Stratified Random Surveys (StRS) of Coral Demography (Adult and Juvenile Corals) across American Samoa in 2015 (NCEI Accession 0159173). NOAA National Centers for Environmental Information. Dataset. https://www.ncei.noaa.gov/archive/accession/0159173.

Coral Reef Ecosystem Program; Pacific Islands Fisheries Science Center (2017). National Coral Reef Monitoring Program: Stratified Random Surveys (StRS) of Coral Demography (Adult and Juvenile Corals) across Jarvis Island from 2016-05-16 to 2016-05-22 (NCEI Accession 0159164). NOAA National Centers for Environmental Information. Dataset. https://www.ncei.noaa.gov/archive/accession/0159164.

Coral Reef Ecosystem Program; Pacific Islands Fisheries Science Center (2017). National Coral Reef Monitoring Program: Stratified Random Surveys (StRS) of Coral Demography (Adult and Juvenile Corals) across the Main Hawaiian Islands from 2013-08-02 to 2013-10-29 (NCEI Accession 0159147). NOAA National Centers for Environmental Information. Dataset. https://accession.nodc.noaa.gov/0159147.

Coral Reef Ecosystem Program; Pacific Islands Fisheries Science Center (2017). National Coral Reef Monitoring Program: Stratified Random Surveys (StRS) of Coral Demography (Adult and Juvenile Corals) across the Mariana Archipelago from 2017-05-05 to 2017-06-21 (NCEI Accession 0166383). NOAA National Centers for Environmental Information. Dataset. https://www.ncei.noaa.gov/archive/accession/0166383.

Ecosystem Sciences Division, Pacific Islands Fisheries Science Center; Papahānaumokuākea Marine National Monument (2019). National Coral Reef Monitoring Program: Stratified random surveys (StRS) of coral demography (adult and juvenile corals) at French Frigate Shoals, Lisianski Island, and Midway Atoll of the Northwestern Hawaiian Islands from 2014-08-14 to 2014-08-26 (NCEI Accession 0184908). NOAA National Centers for Environmental Information. Dataset. https://www.ncei.noaa.gov/archive/accession/0184908. 

Coral Reef Ecosystem Program; Pacific Islands Fisheries Science Center (2017). National Coral Reef Monitoring Program: Stratified Random Surveys (StRS) of Coral Demography (Adult and Juvenile Corals) across the Pacific Remote Island Areas from 2015-01-26 to 2015-04-28 (NCEI Accession 0159161). NOAA National Centers for Environmental Information. Dataset. https://accession.nodc.noaa.gov/0159161.

Coral Reef Ecosystem Program; Pacific Islands Fisheries Science Center (2017). National Coral Reef Monitoring Program: Stratified Random Surveys (StRS) of Coral Demography (Adult and Juvenile Corals) across Rose Atoll from 2016-05-01 to 2016-05-04 (NCEI Accession 0159167). NOAA National Centers for Environmental Information. Dataset. https://www.ncei.noaa.gov/archive/accession/0159167. 

Coral Reef Ecosystem Program; Pacific Islands Fisheries Science Center (2017). National Coral Reef Monitoring Program: Stratified Random Surveys (StRS) of Coral Demography (Adult and Juvenile Corals) across Wake Island from 2014-03-16 to 2014-03-20 (NCEI Accession 0159162). NOAA National Centers for Environmental Information. Dataset. https://www.ncei.noaa.gov/archive/accession/0159162.