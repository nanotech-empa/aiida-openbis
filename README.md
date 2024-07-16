# D2.1 - Metadata schema to store and access microscopy data
Common metadata schema to store and access microscopy data from simulations and experiments in a single platform.

## Authors
- Aliaksandr Yakutovich
- Fabio Lopes
- Carlo Pignedoli

## Goal
* `Notebooks` contains schemas and Pythron scripts needed to set up openBIS instance.
* `setup_openbis_using_linkml.py` creates vocabularies, property types, and object types following the schema provided in the `materialMLinfo.yaml` file contained in the `Metadata_Schemas_LinkML` directory. It also creates the spaces, projects, and collections that will contain data inserted by users. To build such elements, it follows the configuration avaialble in the `collection_config.json`.
* `materialMLinfo.yaml` is a YAML file following the [LinkML](https://linkml.io/) file format. It contains the schema required to mimic what is done in the material science labs working with molecules on surfaces and scanning probe microscopy. It also contains the ontology links necessary to make the data FAIR.
* `import_crystals_to_openBIS.ipynb` uploads the information stored in the inventory of crystals into openBIS.
* `import_molecules_to_openBIS.ipynb` uploads the information stored in the inventory of molecules into openBIS.
* `import_layered2dmaterials_to_openBIS.ipynb` uploads the information stored in the inventory of 2D-layer materials into openBIS.
* `import_images_to_openbis.ipynb` uploads SXM and DAT files into openBIS imaging plugin.

## Achievement
This repository contains all the files needed for setting up the openBIS instance.

## External links
- https://www.ebi.ac.uk/ols4/
- https://linkml.io/

## Acknowledgements
The [PREMISE](https://ord-premise.github.io/) project is supported by the [Open Research Data Program](https://ethrat.ch/en/eth-domain/open-research-data/) of the ETH Board.

![image](https://github.com/ord-premise/metadata-batteries/assets/45081142/74640b5c-ee94-41e1-9acd-fa47da866fe8)