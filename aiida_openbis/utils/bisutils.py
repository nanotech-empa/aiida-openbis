from pybis import Openbis
from aiida.orm import load_node
import tempfile
import shutil
from aiida.common import NotExistent

from rdkit import Chem  # pylint: disable=(import-error)
from rdkit.Chem import AllChem  # pylint: disable=(import-error)


def log_in(
    bisurl='openbis-empa-lab205.labnotebook.ch/openbis/webapp/eln-lims/',
    bisuser='cpignedoli',
    bispasswd='openbisempa'
):
    """Function to login to openBIS."""
    session = Openbis(bisurl)
    session.login(bisuser, bispasswd, save_token=True)
    return session


def log_out(session=None):
    """Function to logout from openBIS."""
    if session and session.is_session_active():
        session.logout()
    return not session.is_session_active()


def allspaces(session=None):
    """Function to retrieve from openBIS all spaces."""
    result = []
    if session and session.is_session_active():
        result = [each.code for each in session.get_spaces()]

    return result


def allprojects(session=None):
    """Function to retrieve from openBIS all projects."""
    result = []
    if session and session.is_session_active():
        result = [each.identifier for each in session.get_projects()]

    return result


def allexperiments(session=None):
    """Function to retrieve from openBIS all experiments."""
    if session and session.is_session_active():
        result = [each.identifier for each in session.get_experiments()]

    return result


def get_molecules(session=None):
    """Function to retrieve from openBIS objects in Molecules collection."""
    result = []
    if session and session.is_session_active():
        available_molecules = session.get_collection('/MATERIALS/SAMPLES/MOLECULES').get_samples()
        result = [(mol.props['$name'], mol.props['molecule.id'], mol.props['molecule.smile'])
                  for mol in available_molecules]
    return result


def new_optimized_geo(session=None, structure=None):
    """Function to export to openBIS an optimized geometry."""
    obj = session.new_object(collection='/MATERIALS/SAMPLES/COMPUTEDGEO', type='STRUCTUREDATA')
    obj.props['structuredata.info'] = 'Optimized geometry'
    obj.save()
    if structure:
        tmpdir = tempfile.mkdtemp()
        file_path1 = tmpdir + "/" + 'struc.xyz'
        file_path2 = tmpdir + "/" + 'struc.png'
        structure.write(file_path1)
        structure.write(file_path2)
        rawds = session.new_dataset(type='RAW_DATA', object=obj, file=file_path1)
        rawds.save()

        preview = session.new_dataset(type='ELN_PREVIEW', object=obj, file=file_path2)
        preview.save()
        shutil.rmtree(tmpdir)
    return obj


def new_molecule(session=None, name=None, molid=None, smile=None, attachment=None):
    """Function  to create in openBIS a new MOLECULE object."""
    obj = session.new_object(collection='/MATERIALS/SAMPLES/MOLECULES', type='MOLECULE')
    obj.props['$name'] = name
    obj.props['molecule.smile'] = smile
    if molid:
        obj.props['molecule.id'] = molid
    obj.save()
    if attachment:
        rawds = session.new_dataset(type='RAW_DATA', object=obj, file=attachment)
        rawds.save()
    tmpl = Chem.MolFromSmiles(smile)
    AllChem.Compute2DCoords(tmpl)
    img = Chem.Draw.MolToImage(tmpl)
    tmpdir = tempfile.mkdtemp()
    file_path = tmpdir + "/" + 'struc.png'
    img.save(file_path)
    preview = session.new_dataset(type='ELN_PREVIEW', object=obj, file=file_path)
    preview.save()

    return obj


def new_product(session=None, name=None, smile=None, theyield=None, length=None, temperature=None):  # pylint: disable=(too-many-arguments)
    """Function  to create in openBIS a new MOLPRODUCT object."""
    obj = session.new_object(collection='/MATERIALS/SAMPLES/PRODUCTS', type='MOLPRODUCT')
    obj.props['$name'] = name
    obj.props['molproduct.smile'] = smile
    if theyield:
        obj.props['molproduct.yield'] = theyield
    if length:
        obj.props['molproduct.length'] = length
    if temperature:
        obj.props['molproduct.temperature'] = temperature
    obj.save()
    tmpl = Chem.MolFromSmiles(smile)
    AllChem.Compute2DCoords(tmpl)
    img = Chem.Draw.MolToImage(tmpl)
    tmpdir = tempfile.mkdtemp()
    file_path = tmpdir + "/" + 'struc.png'
    img.save(file_path)
    preview = session.new_dataset(type='ELN_PREVIEW', object=obj, file=file_path)
    preview.save()
    return obj


def new_reaction_products(reactions=None, molecules=None, attachment=None):
    """Function  to create in openBIS objects from a cdxml reaction."""
    session = log_in()
    allm = set(m['name'] for m in molecules) - set(p['product'] for p in reactions)
    allobj = {}
    for mol in molecules:
        if mol['name'] in allm:
            allobj[mol['name']] = new_molecule(
                session=session, name=mol['name'], smile=mol['smile'], attachment=attachment
            ).permId

        else:
            allobj[mol['name']] = new_product(session=session, name=mol['name'], smile=mol['smile']).permId

    for reac in reactions:
        print('adding ', reac['product'], ' id ', allobj[reac['product']])
        print('as child of ', reac['reactant'], ' id ', allobj[reac['reactant']])
        #print('parent is new ', allobj[reac['reactant']].is_new)
        #print('parent is new ', allobj[reac['product']].is_new)
        reactant = session.get_object(allobj[reac['reactant']])
        product = session.get_object(allobj[reac['product']])
        #allobj[reac['reactant']].save()
        #allobj[reac['reactant']].add_children(allobj[reac['product']])
        #allobj[reac['reactant']].save()
        if reac['yield']:
            product.props['molproduct.yield'] = reac['yield']
        if reac['length']:
            product.props['molproduct.length'] = reac['length']
        if reac['temperature']:
            product.props['molproduct.temperature'] = reac['temperature']
        product.save()
        reactant.add_children(product)
        reactant.save()
    session.logout()
    return not session.is_session_active()


def aiidalab_geo_opt(pk=None, collection='/SPIN_CHAIN/TRIANGULENE_BASED/TRIANGULENE_BASED_EXP_2'):
    """Function to export to openBIS results from an AiiDA geo opt workflow."""
    if pk:
        try:
            node = load_node(pk)
        except NotExistent:
            return False
    # Open openBIS session.
    session = log_in()
    # Create the new sstructure data openBIS object.
    asegeo = node.outputs.output_structure.get_ase()

    newSD = new_optimized_geo(session=session, structure=asegeo)
    newgeoopt = session.new_object(collection=collection, type='AIIDALAB_GEO_OPT')
    newgeoopt.add_children(newSD)
    newgeoopt.add_parents('/MATERIALS/ORGANIZATION/PER1')
    newgeoopt.props['start_date'] = node.ctime.strftime("%Y-%m-%d %H:%M:%S")
    newgeoopt.props['end_date'] = node.mtime.strftime("%Y-%m-%d %H:%M:%S")
    newgeoopt.props['aiidalab_geo_opt.description'] = node.description
    newgeoopt.props['aiidalab_geo_opt.pk'] = node.pk
    newgeoopt.props['aiidalab_geo_opt.uuid'] = node.uuid
    newgeoopt.props[
        'aiidalab_geo_opt.url'
    ] = 'https://aiidalab.materialscloud.org/user/carlo.pignedoli@empa.ch/apps/apps/home/start.ipynb?'
    newgeoopt.props['aiidalab_geo_opt.outdict'] = str(dict(node.outputs.output_parameters))
    newgeoopt.save()

    # Close openBIS session
    session.logout()
    return True
