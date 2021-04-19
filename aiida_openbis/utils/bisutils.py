from pybis import Openbis
from aiida.orm import load_node
from ase.io import write
import tempfile
import shutil


def log_in(
        bisurl='openbis-empa-lab205.labnotebook.ch/openbis/webapp/eln-lims/',
        bisuser='cpignedoli',
        bispasswd='openbisempa'):
    session = Openbis(bisurl)
    session.login(bisuser, bispasswd, save_token=True)
    return session


def log_out(session=None):
    if session and session.is_session_active():
        session.logout()
    return


def allspaces(session=None):
    if session and session.is_session_active():
        return [each.code for each in session.get_spaces()]


def allprojects(session=None):
    if session and session.is_session_active():
        return [each.identifier for each in session.get_projects()]


def allexperiments(session=None):
    if session and session.is_session_active():
        return [each.identifier for each in session.get_experiments()]


def get_molecules(session=None):
    if session and session.is_session_active():
        available_molecules = session.get_collection(
            '/MATERIALS/SAMPLES/MOLECULES').get_samples()
        return [(mol.props['$name'], mol.props['molecule.id'],
                 mol.props['molecule.smile']) for mol in available_molecules]


def new_optimized_geo(session=None, structure=None):
    obj = session.new_object(collection='/MATERIALS/SAMPLES/COMPUTEDGEO',
                             type='STRUCTUREDATA')
    obj.props['structuredata.info'] = 'Optimized geometry'
    obj.save()
    if structure:
        tmpdir = tempfile.mkdtemp()
        file_path1 = tmpdir + "/" + 'struc.xyz'
        file_path2 = tmpdir + "/" + 'struc.png'
        structure.write(file_path1)
        structure.write(file_path2)
        rawds = session.new_dataset(type='RAW_DATA',
                                    object=obj,
                                    file=file_path1)
        rawds.save()

        preview = session.new_dataset(type='ELN_PREVIEW',
                                      object=obj,
                                      file=file_path2)
        preview.save()
        shutil.rmtree(tmpdir)
    return obj


def new_molecule(session=None,
                 name=None,
                 molid=None,
                 smile=None,
                 attachment=None):
    obj = session.new_object(collection='/MATERIALS/SAMPLES/MOLECULES',
                             type='MOLECULE')
    obj.props['$name'] = name
    obj.props['molecule.smile'] = smile
    obj.props['molecule.id'] = molid
    obj.save()
    if attachment:
        rawds = session.new_dataset(type='RAW_DATA',
                                    object=obj,
                                    file=attachment)
        rawds.save()
    return obj


def new_product(session=None, name=None, molid=None, smile=None):
    obj = session.new_object(collection='/MATERIALS/SAMPLES/PRODUCTS',
                             type='MOLPRODUCT')
    obj.props['$name'] = name
    obj.props['molecule.smile'] = smile
    obj.props['molecule.id'] = molid
    obj.save()
    return obj


def new_reaction_products(session=None,
                          reactions=None,
                          molecules=None,
                          attachment=None):
    session = log_in()
    allm = set([m['name']
                for m in molecules]) - set([p['product'] for p in reactions])
    allobj = {}
    for mol in molecules:
        if mol['name'] in allm:
            allobj[mol['name']] = new_molecule(session=session,
                                               name=mol['name'],
                                               smile=mol['smile'],
                                               attachment=attachment)

        else:
            allobj[mol['name']] = new_product(session=session,
                                              name=mol['name'],
                                              smile=mol['smile'])

    for reac in reactions:
        allobj[reac['reactant']].add_children(allobj[reac['product']])
        allobj[reac['reactant']].save()
    session.logout()


def aiidalab_geo_opt(
        pk=None,
        collection='/SPIN_CHAIN/TRIANGULENE_BASED/TRIANGULENE_BASED_EXP_2'):
    if pk:
        try:
            node = load_node(pk)
        except:
            return
        # Open openBIS session.
        session = log_in()
        # Create the new sstructure data openBIS object.
        asegeo = node.outputs.output_structure.get_ase()

        newSD = new_optimized_geo(session=session, structure=asegeo)
        newgeoopt = session.new_object(collection=collection,
                                       type='AIIDALAB_GEO_OPT')
        newgeoopt.add_children(newSD)
        newgeoopt.add_parents('/MATERIALS/ORGANIZATION/PER1')
        newgeoopt.props['start_date'] = node.ctime.strftime(
            "%Y-%m-%d %H:%M:%S")
        newgeoopt.props['end_date'] = node.mtime.strftime("%Y-%m-%d %H:%M:%S")
        newgeoopt.props['aiidalab_geo_opt.description'] = node.description
        newgeoopt.props['aiidalab_geo_opt.pk'] = node.pk
        newgeoopt.props['aiidalab_geo_opt.uuid'] = node.uuid
        newgeoopt.props[
            'aiidalab_geo_opt.url'] = 'https://aiidalab.materialscloud.org/user/carlo.pignedoli@empa.ch/apps/apps/home/start.ipynb?'
        newgeoopt.props['aiidalab_geo_opt.outdict'] = str(
            dict(node.outputs.output_parameters))
        newgeoopt.save()

        # Close openBIS session
        session.logout()
