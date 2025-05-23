# Sporophere_Bacteria

Influence of AMF spore traits on associated bacteria

___


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

### Taxonomic database

Silva nr99 database v [138.2](https://zenodo.org/records/14169026)

Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO (2013) The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucl. Acids Res. 41 (D1): D590-D596.

___





1. Match up ./data/SeqCoast_Manifest.xlsx with sample metadata from Bala's Google Drive for this project
2. Build single complete metadata sheet
3. Run all samples through dada2 pipeline (use silva taxonomy?)
4. Remove contaminants (identified using negative controls)

5. background soil *should* resemble the spore surface bacteria for each AMF taxon
6. 

___


