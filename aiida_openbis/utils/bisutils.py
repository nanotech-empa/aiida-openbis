from pybis import Openbis
from aiida.orm import load_node
import tempfile
import zipfile
import shutil
from aiida.common import NotExistent

from rdkit import Chem  # pylint: disable=(import-error)
from rdkit.Chem import AllChem  # pylint: disable=(import-error)


def log_in(bisurl='openbis', bisuser='admin', bispasswd='changeit'):
    """Function to login to openBIS."""
    if Openbis(bisurl, verify_certificates=False).is_token_valid():
        session = Openbis(bisurl, verify_certificates=False)
    else:
        Openbis(bisurl, verify_certificates=False).login(bisuser, bispasswd, save_token=True)
        session = Openbis(bisurl, verify_certificates=False)
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
        result = [
            (mol.props['name'], mol.permId, mol.props['smiles']) for mol in available_molecules
        ]
    return result


def get_precursors(session=None):
    """Function to retrieve from openBIS objects in Molecules collection."""
    result = []
    if session and session.is_session_active():
        available_molecules = session.get_collection('/MATERIALS/SAMPLES/MOLECULES').get_samples()
        result = [(mol.permId, mol.props['name'], mol.props['id'], mol.props['smiles'])
                  for mol in available_molecules]
    return result


def new_chem_sketch(session=None, attachment=None):
    if not attachment:
        return False
    obj = session.new_object(collection='/MATERIALS/SAMPLES/CHEMSKETCH', type='CHEMSKETCH')
    obj.props['description'] = 'Chemical sketch of ....'
    obj.save()
    if attachment:
        rawds = session.new_dataset(type='RAW_DATA', object=obj, file=attachment)
        rawds.save()
    return obj


def new_optimized_geo(session=None, pk=None):
    """Function to export to openBIS an optimized geometry."""
    node = load_node(pk)
    structure = node.get_ase()
    obj = session.new_object(collection='/MATERIALS/SAMPLES/COMPUTEDGEO', type='STRUCTUREDATA')
    obj.props['description'] = 'Optimized geometry'
    obj.props['aiida.uuid'] = node.uuid
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


def get_opt_geo_ids(session=None):
    return [(obj.permId, obj.props['aiida.uuid'])
            for obj in session.get_objects(collection='/MATERIALS/SAMPLES/COMPUTEDGEO')]


def new_molecule(session=None, name=None, molid=None, smile=None, cdxml=None):
    """Function  to create in openBIS a new MOLECULE object."""
    obj = session.new_object(collection='/MATERIALS/SAMPLES/MOLECULES', type='MOLECULE')
    obj.props['name'] = name
    obj.props['smiles'] = smile
    if molid:
        obj.props['id'] = molid
    obj.save()
    tmpl = Chem.MolFromSmiles(smile)
    AllChem.Compute2DCoords(tmpl)
    img = Chem.Draw.MolToImage(tmpl)
    tmpdir = tempfile.mkdtemp()
    file_path = tmpdir + "/" + 'struc.png'
    img.save(file_path)
    preview = session.new_dataset(type='ELN_PREVIEW', object=obj, file=file_path)
    preview.save()
    if cdxml:
        tmpdir = tempfile.mkdtemp()
        file_path = tmpdir + "/" + 'sketch.cdxml'
        with open(file_path, 'w') as f:
            f.write(cdxml)
        rawds = session.new_dataset(type='RAW_DATA', object=obj, file=file_path)
        rawds.save()

    return obj


def new_molecule_precursor(
    session=None,
    number=None,
    batch=None,
    acronym=None,
    iupac=None,
    chemformula=None,
    smile=None,
    cas=None,
    evaporationt=None,
    hazardous=None,
    azardousdescription=None,
    fridge=None,
    nolight=None,
    dry=None,
    nooxygen=None,
    otherstorage=None,
    specifyotherstorage=None,
    supplier=None,
    synthesizedby=None,
    supplierownname=None,
    amount=None,
    receivingdate=None,
    cdxml=None,
    addcomments=None,
    filestoupload=None
):
    """Function  to create in openBIS a new MOLECULE object."""
    session = log_in()
    obj = session.new_object(collection='/MATERIALS/MOLECULES/PRECURSORS', type='MOLECULE')
    obj.props['molecule.number'] = number
    obj.props['molecule.batch'] = batch
    obj.props['molecule.acronym'] = acronym
    obj.props['molecule.iupac'] = iupac
    obj.props['molecule.chemformula'] = chemformula
    obj.props['molecule.smile'] = smile
    obj.props['molecule.cas'] = cas
    obj.props['molecule.evaporationt'] = evaporationt
    obj.props['molecule.hazardous'] = hazardous
    obj.props['molecule.azardousdescription'] = azardousdescription
    obj.props['molecule.fridge'] = fridge
    obj.props['molecule.nolight'] = nolight
    obj.props['molecule.dry'] = dry
    obj.props['molecule.nooxygen'] = nooxygen
    obj.props['molecule.otherstorage'] = otherstorage
    obj.props['molecule.specifyotherstorage'] = specifyotherstorage
    obj.props['molecule.supplier'] = supplier
    obj.props['molecule.synthesizedby'] = synthesizedby
    obj.props['molecule.supplierownname'] = supplierownname
    obj.props['molecule.amount'] = str(amount)
    obj.props['molecule.receivingdate'] = receivingdate
    obj.props['molecule.addcomments'] = addcomments

    obj.save()
    tmpl = Chem.MolFromSmiles(smile)
    AllChem.Compute2DCoords(tmpl)
    img = Chem.Draw.MolToImage(tmpl)
    tmpdir = tempfile.mkdtemp()
    file_path = tmpdir + "/" + 'struc.png'
    img.save(file_path)
    print('uploading molecule preview')
    preview = session.new_dataset(type='ELN_PREVIEW', object=obj, file=file_path)
    preview.save()
    if cdxml:
        tmpdir = tempfile.mkdtemp()
        file_path = tmpdir + "/" + 'sketch.cdxml'
        with open(file_path, 'w') as f:
            f.write(cdxml)
        rawds = session.new_dataset(type='RAW_DATA', object=obj, file=file_path)
        rawds.save()
    for addfile in filestoupload:
        print('uploading ', addfile)
        rawds = session.new_dataset(type='RAW_DATA', object=obj, file=addfile)
        rawds.save()

    return obj


def new_product(
    session=None,
    name=None,
    smile=None,
    cdxml=None,
    theyield=None,
    length=None,
    temperature=None
):  # pylint: disable=(too-many-arguments)
    """Function  to create in openBIS a new PRODUCT object."""
    obj = session.new_object(collection='/MATERIALS/SAMPLES/PRODUCTS', type='MOLECULE_PRODUCT')
    obj.props['name'] = name
    obj.props['smiles'] = smile
    if theyield:
        obj.props['yield'] = theyield
    if length:
        obj.props['length'] = length
    if temperature:
        obj.props['temperature'] = temperature
    obj.save()
    tmpl = Chem.MolFromSmiles(smile)
    AllChem.Compute2DCoords(tmpl)
    img = Chem.Draw.MolToImage(tmpl)
    tmpdir = tempfile.mkdtemp()
    file_path = tmpdir + "/" + 'struc.png'
    img.save(file_path)
    preview = session.new_dataset(type='ELN_PREVIEW', object=obj, file=file_path)
    preview.save()
    if cdxml:
        tmpdir = tempfile.mkdtemp()
        file_path = tmpdir + "/" + 'sketch.cdxml'
        with open(file_path, 'w') as f:
            f.write(cdxml)
        rawds = session.new_dataset(type='RAW_DATA', object=obj, file=file_path)
        rawds.save()
    return obj


def new_reaction_products(reactions=None, molecules=None, attachment=None):
    """Function  to create in openBIS objects from a cdxml reaction."""
    if not attachment:
        return False

    session = log_in()
    cdxmlid = new_chem_sketch(session=session, attachment=attachment).permId
    allm = set(m['name'] for m in molecules) - set(p['product'] for p in reactions)
    allobj = {}
    for mol in molecules:
        if mol['name'] in allm:
            allobj[mol['name']] = new_molecule(
                session=session, name=mol['name'], smile=mol['smile'], cdxml=mol['cdxml']
            ).permId

        else:
            allobj[mol['name']] = new_product(
                session=session, name=mol['name'], smile=mol['smile'], cdxml=mol['cdxml']
            ).permId

    for reac in reactions:
        reactant = session.get_object(allobj[reac['reactant']])
        product = session.get_object(allobj[reac['product']])
        cdxml = session.get_object(cdxmlid)
        if reac['yield']:
            product.props['yield'] = reac['yield']
        if reac['length']:
            product.props['length'] = reac['length']
        if reac['temperature']:
            product.props['temperature'] = reac['temperature']
        product.save()
        cdxml.add_children(product)
        cdxml.save()
        reactant.add_children(product)
        reactant.save()
        # The reaction .cdxml is children of precursor molecules.
        if reac['reactant'] in allm:
            reactant.add_children(cdxml)
            reactant.save()
        # Product molecules are children of the reaction .cdxml
        else:
            cdxml.add_children(reactant)
            cdxml.save()
    session.logout()
    return not session.is_session_active()


def is_from_openbis(pk=None):
    permId = False
    if pk:
        try:
            node = load_node(pk)
        except NotExistent:
            return False
        try:
            permId = node.extras['eln']['sample_uuid']
        except KeyError:
            permId = False
    return permId


def aiidalab_geo_opt(session=None, pk=None, collection='/DEFAULT_LAB_NOTEBOOK/7-AGNR/PRISTINE'):
    """Function to export to openBIS results from an AiiDA geo opt workflow."""
    if pk:
        try:
            node = load_node(pk)
        except NotExistent:
            return False
    # Open openBIS session.
    close_at_end = False
    if session is None:
        session = log_in()
        close_at_end = True
    # Create the new sstructure data openBIS object.
    structure_pk = node.outputs.output_structure.pk
    # Check if the input workchain input structure originated from openBIS and gets permId
    bis_parent_permId = is_from_openbis(
        pk=node.inputs.cp2k__file__input_xyz.creator.inputs.structure.pk
    )
    newSD = new_optimized_geo(session=session, pk=structure_pk)
    newgeoopt = session.new_object(collection=collection, type='AIIDALAB_GEO_OPT')
    newgeoopt.add_children(newSD)
    #newgeoopt.add_parents('/MATERIALS/ORGANIZATION/PER1')
    newgeoopt.props['start_date'] = node.ctime.strftime("%Y-%m-%d %H:%M:%S")
    newgeoopt.props['end_date'] = node.mtime.strftime("%Y-%m-%d %H:%M:%S")
    newgeoopt.props['description'] = node.description
    newgeoopt.props['aiida.uuid'] = node.uuid
    newgeoopt.props['aiida.server_url'] = 'http://127.0.0.1:8888/apps/apps/home/start.ipynb'
    newgeoopt.save()
    if bis_parent_permId:
        newgeoopt.add_parents(bis_parent_permId)
        newgeoopt.save()

    # Close openBIS session
    info = [newSD.permId, newgeoopt.permId]
    if close_at_end:
        session.logout()
    return info


def aiidalab_spm(zip_path=None, pk=None, collection='/DEFAULT_LAB_NOTEBOOK/7-AGNR/PRISTINE'):  # pylint: disable=(too-many-locals)
    """Function to export to openBIS STM sets from an AiiDAlab SPM workflow."""
    if pk:
        try:
            node = load_node(pk)
        except NotExistent:
            return False
    # Open openBIS session.
    session = log_in()
    newspm = session.new_object(collection=collection, type='AIIDALAB_SPM')

    #newspm.add_parents(the_geo)
    # identify user and add as parent
    #newspm.add_parents('/MATERIALS/ORGANIZATION/PER1')
    newspm.props['start_date'] = node.ctime.strftime("%Y-%m-%d %H:%M:%S")
    newspm.props['end_date'] = node.mtime.strftime("%Y-%m-%d %H:%M:%S")
    newspm.props['description'] = node.description
    newspm.props['aiida.uuid'] = node.uuid
    newspm.props['aiida.server_url'] = 'https://aiidalab.materialscloud.org/user/' + node.user.email
    newspm.props['aiida.notes'] = 'Notes:'
    newspm.save()
    # search parent geometry and add as parent
    # if not present check if originated from geo_opt and add geo_opt
    structure_uuid = node.inputs.structure.uuid
    available_geos = get_opt_geo_ids(session=session)
    try:
        idn = [a[1] for a in available_geos].index(structure_uuid)
        structure_permId = available_geos[idn][0]
    except ValueError:
        try:
            structure_permId = aiidalab_geo_opt(
                session=session, pk=node.inputs.structure.creator.caller.pk
            )[0]
        except (AttributeError, ValueError):
            structure_permId = None
    if structure_permId:
        newspm.add_parents(structure_permId)
        newspm.save()
    # Attach zipfile.
    rawds = session.new_dataset(type='RAW_DATA', object=newspm, file=zip_path)
    rawds.props['$name'] = 'Igor_files'
    rawds.props['notes'] = 'Zip file containing raw data in .igor and .txt format'
    rawds.save()

    # Gallery of STM images
    xmlstring = '<?xml version="1.0" encoding="UTF-8"?>\n<html><head>Images</head><body>'
    thezip = zipfile.ZipFile(zip_path, 'r')
    # Parse through the files.
    for filename in thezip.namelist():
        if filename.endswith('.png'):
            pngfile = tempfile.mkdtemp() + "/" + filename
            # save the image as raw_dataset
            with open(pngfile, 'wb') as newf:
                newf.write(thezip.open(filename).read())
                rawds = session.new_dataset(type='RAW_DATA', object=newspm, file=pngfile)
                rawds.props['$name'] = filename.replace('.png', '')
                rawds.props['notes'] = 'SPM png file'
                rawds.save()
            img_url = 'https://localhost:8443/datastore_server/'
            img_url += str(rawds.permId) + '/' + list(rawds.file_links.keys())[0]
            # New image.
            xmlstring += '<figure class="image image-style-align-left image_resized" style="width:50.0%;">'
            xmlstring += '<img src="' + img_url + '">'
            # Figure caption.
            xmlstring += '<figcaption>Parameters: '
            xmlstring += filename.replace('.png', '')
            xmlstring += '</figcaption></figure>'
    xmlstring += '</body></html>'
    newspm.props['xml.image'] = xmlstring
    newspm.save()
    # Close openBIS session
    session.logout()
    return True
