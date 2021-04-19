from pybis import Openbis


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


bisdata = log_in()
#print(aa.is_session_active())
##print(allspaces(session=aa))
#print(allprojects(session=aa))
#print(allexperiments(session=aa))
print(get_molecules(session=bisdata))
log_out(session=bisdata)
print(bisdata.is_session_active())
#print(aa.is_session_active())
