# from flask import Flask, redirect, url_for, render_template, request, session, flash
# from datetime import timedelta
# from molecule import molecule

# app = Flask(__name__)
# app.register_blueprint(molecule, url_prefix="/molecule")

# @app.route("/")
# def home():
#     return render_template("index.html")

# @app.route('/new_space')
# def new_space():
#     return render_template('new_space.html')


# if __name__ == "__main__":
#     app.run(debug = True)

from fastapi import FastAPI, Request
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates
from fastapi.staticfiles import StaticFiles
from aiida_openbis.utils import bisutils
import re
import uvicorn
from molecule import molecule_router
from new_space import new_space_router
from new_project import new_project_router

# Create a FastAPI instance
app = FastAPI()

# Configure templates for rendering HTML
templates = Jinja2Templates(directory="templates")

# Mount the "static" directory for serving static files
app.mount("/static", StaticFiles(directory="static"), name="static")

app.include_router(molecule_router)
app.include_router(new_space_router)
app.include_router(new_project_router)

# Create a session
session = bisutils.log_in(bisurl="openbis", bisuser="admin", bispasswd="changeit")

@app.get("/", response_class=HTMLResponse)
async def home(request: Request):
    return templates.TemplateResponse("index.html", {"request": request})

@app.get("/new_space", response_class=HTMLResponse)
async def new_space(request: Request):
    return templates.TemplateResponse("new_space.html", {"request": request})

@app.get("/new_project", response_class=HTMLResponse)
async def new_project(request: Request):
    all_spaces = session.get_spaces().df.code
    return templates.TemplateResponse("new_project.html", {"request": request, "items": all_spaces})

@app.get("/new_experiment", response_class=HTMLResponse)
async def new_experiment(request: Request):
    all_spaces = session.get_spaces().df.code
    return templates.TemplateResponse("new_experiment.html", {"request": request, "items": all_spaces})

@app.get("/new_collection", response_class=HTMLResponse)
async def new_collection(request: Request):
    all_projects = session.get_projects().df.identifier
    all_objects_types = session.get_object_types().df.code
    return templates.TemplateResponse("new_collection.html", {"request": request, "projects_codes": all_projects,
                                                              "objects_types": all_objects_types})


if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=5000)