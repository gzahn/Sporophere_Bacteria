# Sporophere_Bacteria

Influence of AMF spore traits on associated bacteria

___

![](https://raw.githubusercontent.com/gzahn/Sporophere_Bacteria/refs/heads/main/output/figs/spore_and_bacteria.png)

## Experimental overview

There are 6 replicate samples of spores isolated from 28 INVAM cultured AMF isolates.

Each pool of spores was divided into 6 piles in a watch glass under a dissecting scope.

For a given isolate (_e.g._, MG101B-12), there are two treatment conditions:

        1. sterile (surface sterilized with 10% bleach and 70% ethanol)
        2. unsterile (not surface sterilized)

3 of the piles were surface sterilized, and 3 were not.

The *total* sporosphere bacterial community can be derived from the *unsterile* replicates of a given isolate. However, this may include 
bacteria present *inside* of spores. The surface sterilized replicates should have bacteria that was *only* inside of spores.
Therefore, the surface bacterial community might be inferred as:

		surface_community = total_community - internal_community
				(unsterilized)  -  (sterlized)

Additionally, bulk soil from the INVAM cultures was sequenced and are identified as treatment == "invam bulk soil".

These can provide the background bacterial community present in the soil.
Therefore the bacterial community that is especially associated with spore surfaces can be inferred as:

                surface_community = total_community - internal_community - invam_soil

Finally, we sequenced 54 soil samples from NEON sites with a range of *Acaulospora* abundances. The idea is that if we see any bacterial community signatures in spores with high surface ornamentation, we might be able to look back at the NEON soil samples and look for that same community signature across a gradient of "ornamentation abundance."

___


### Taxonomic database

Silva nr99 database v [138.2](https://zenodo.org/records/14169026)

Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO (2013) The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucl. Acids Res. 41 (D1): D590-D596.

___





- Run all samples through dada2 pipeline (use silva taxonomy?)
- Remove contaminants (identified using negative controls)
- Null: background soil *should* resemble the spore surface bacteria for each AMF taxon


___


