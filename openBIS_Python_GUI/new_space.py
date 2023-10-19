from fastapi import APIRouter, Request
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.templating import Jinja2Templates
from aiida_openbis.utils import bisutils
import re
import utils

new_space_router = APIRouter()

templates = Jinja2Templates(directory="templates")

session = bisutils.log_in(bisurl="openbis", bisuser="admin", bispasswd="changeit")

# Define a FastAPI route that corresponds to the "/molecule" blueprint in Flask
@new_space_router.post("/create_new_space")
async def create_new_space(data: dict):
    space_data = data.get('data')
    space_name = space_data.get('spaceName')
    space_description = space_data.get('spaceDescription')

    # Convert the name to a code that openBIS API can understand
    new_space_code = utils.convert_name_to_code(space_name)

    # Check if the space it is being created already exists in the DB
    new_space_not_exists = session.get_spaces(code = new_space_code).df.empty

    response = {}
    if new_space_not_exists:
        values = {'description': space_description}
        session.new_space(code=new_space_code, **values).save()
        response['result'] = f'Space {space_name} created successfully!'
    else:
        response['result'] = f'Space {space_name} already exists!'

    return JSONResponse(content=response)