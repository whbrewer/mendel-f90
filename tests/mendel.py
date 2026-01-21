#!/usr/bin/env python
for i in range(ng):
    gp.init()
    # randomly select individuals and send their genetic information to the 
    # appropriate neighboring tribes
	#gp.migrate()
	for j in range(np/2):
        # randomly mate half of the population with members from other half
        # offspring receives half its genetic makeup from each of its two parents; 
		# add new random mutations to offspring genome
		gp.mate()
	# impose selection based on phenotypic fitness to reduce the population size
	gp.select()
