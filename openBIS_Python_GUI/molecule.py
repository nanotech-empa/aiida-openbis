# from flask import Blueprint, render_template, request, jsonify
# from aiida_openbis.utils import bisutils
# import re

# session = bisutils.log_in(bisurl="openbis", bisuser="admin", bispasswd="changeit")
# molecule = Blueprint("molecule", __name__, static_folder = "static", template_folder = "templates")

# @molecule.route("/")
# def load_molecule_page():
#     return render_template("molecule.html")

# @molecule.route("/update_molecule", methods = ["POST"])
# def update_molecule():
#     all_molecules = session.get_samples(collection="MOLECULES", props=['empa-number']).df
#     selected_empa_number = request.get_json()
#     selected_empa_number = selected_empa_number.get('data')
#     permId = all_molecules[all_molecules['EMPA-NUMBER']==selected_empa_number]['permId'].values[0]
#     molecule = session.get_sample(permId)

#     response = {}
#     response['moleculeId'] = molecule.code
#     response['moleculeEmpaNumber'] = molecule.props['empa-number']
#     response['moleculeBatch'] = molecule.props['batch']
#     response['moleculeAcronym'] = molecule.props['acronym']
#     response['moleculeIupacName'] = molecule.props['iupac-name']
#     response['moleculeFormula'] = molecule.props['formula']
#     response['moleculeSmiles'] = molecule.props['smiles']
#     response['moleculeCasNumber'] = molecule.props['cas-number']
#     response['moleculeHazardous'] = bool(molecule.props['hazardous'])
#     response['moleculeHazardousSpecify'] = molecule.props.get('hazardous-specify',"")
#     response['moleculeLight'] = bool(molecule.props['light'])
#     response['moleculeOxygen'] = bool(molecule.props['oxygen'])
#     response['moleculeFreezer'] = bool(molecule.props['freezer'])
#     response['moleculeDry'] = bool(molecule.props['dry'])
#     response['moleculeStorageConditionOther'] = molecule.props['storage-condition-other']
#     response['moleculeStorageConditionOtherSpecify'] = molecule.props.get('storage-condition-other-specify',"")
#     response['moleculeAmount'] = molecule.props['amount']
#     response['moleculeReceivingDate'] = re.findall(r'(\d{4}-\d{2}-\d{2})', molecule.props['receiving-date'])[0]
#     response['moleculeComments'] = molecule.props.get('comments', "")
    
#     print(response)
#     return jsonify(response)

from fastapi import APIRouter, Request
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.templating import Jinja2Templates
from aiida_openbis.utils import bisutils
import re

molecule_router = APIRouter()

templates = Jinja2Templates(directory="templates")

session = bisutils.log_in(bisurl="openbis", bisuser="admin", bispasswd="changeit")

# Define a FastAPI route that corresponds to the "/molecule" blueprint in Flask
@molecule_router.get("/molecule", response_class=HTMLResponse)
async def load_molecule_page(request: Request):
    batch_codes = session.get_vocabulary('BATCH-ID').get_terms().df.code
    batch_chars = session.get_vocabulary('BATCH-ID').get_terms().df.label
    return templates.TemplateResponse("molecule.html", {"request": request, "batch_chars":batch_chars, "batch_codes":batch_codes})

@molecule_router.post("/molecule/update_molecule")
async def update_molecule(data: dict):
    all_molecules = session.get_samples(collection="MOLECULES", props=['empa-number']).df
    selected_empa_number = data.get('data')
    permId = all_molecules[all_molecules['EMPA-NUMBER'] == selected_empa_number]['permId'].values[0]
    molecule = session.get_sample(permId)

    response = {
        'moleculeId': molecule.code,
        'moleculeEmpaNumber': molecule.props['empa-number'],
        'moleculeBatch': molecule.props['batch'],
        'moleculeAcronym': molecule.props['acronym'],
        'moleculeIupacName': molecule.props['iupac-name'],
        'moleculeFormula': molecule.props['formula'],
        'moleculeSmiles': molecule.props['smiles'],
        'moleculeCasNumber': molecule.props['cas-number'],
        'moleculeHazardous': bool(molecule.props.get('hazardous')),
        'moleculeHazardousSpecify': molecule.props.get('hazardous-specify', ""),
        'moleculeLight': bool(molecule.props.get('light')),
        'moleculeOxygen': bool(molecule.props.get('oxygen')),
        'moleculeFreezer': bool(molecule.props.get('freezer')),
        'moleculeDry': bool(molecule.props.get('dry')),
        'moleculeStorageConditionOther': molecule.props.get('storage-condition-other', ""),
        'moleculeStorageConditionOtherSpecify': molecule.props.get('storage-condition-other-specify', ""),
        'moleculeAmount': molecule.props['amount'],
        'moleculeReceivingDate': re.findall(r'(\d{4}-\d{2}-\d{2})', molecule.props['receiving-date'])[0],
        'moleculeComments': molecule.props.get('comments', ""),
    }
    
    return JSONResponse(content=response)
