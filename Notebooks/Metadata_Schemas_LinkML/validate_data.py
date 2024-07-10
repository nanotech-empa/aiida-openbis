from linkml.validator import validate

instance = {
    "name": "methane",
    "iupac_name": "methane",
    "sum_formula": "CH4",
    "smiles": "C",
}

report = validate(instance, "materialMLinfo.yaml", "Molecule")

if not report.results:
    print('The instance is valid!')
else:
    for result in report.results:
        print(result.message)