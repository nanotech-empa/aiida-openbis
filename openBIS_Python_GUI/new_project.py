from fastapi import APIRouter, Request
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.templating import Jinja2Templates
from aiida_openbis.utils import bisutils
import re
import utils

new_project_router = APIRouter()

templates = Jinja2Templates(directory="templates")

session = bisutils.log_in(bisurl="openbis", bisuser="admin", bispasswd="changeit")

# Define a FastAPI route that corresponds to the "/molecule" blueprint in Flask
@new_project_router.post("/create_new_project")
async def create_new_project(data: dict):
    project_data = data.get('data')
    project_space = project_data.get('projectSpace')
    project_name = project_data.get('projectName')
    project_description = project_data.get('projectDescription')

    # Convert the name to a code that openBIS API can understand
    new_project_code = utils.convert_name_to_code(project_name)

    # Check if the space it is being created already exists in the DB
    new_space_not_exists = session.get_projects(code = new_project_code, space = project_space).df.empty

    response = {}
    if new_space_not_exists:
        values = {'space': project_space, 'description': project_description}
        session.new_project(code = new_project_code, **values).save()
        response['result'] = f'Project {project_name} created successfully!'
    else:
        response['result'] = f'Project {project_name} already exists!'

    return JSONResponse(content=response)