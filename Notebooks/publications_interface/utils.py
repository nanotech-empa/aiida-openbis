import ipywidgets as ipw
import ipysheet

def Text(**kwargs):
    return ipw.Text(**kwargs)

def Textarea(**kwargs):
    return ipw.Textarea(**kwargs)

def IntText(**kwargs):
    return ipw.IntText(**kwargs)

def FloatText(**kwargs):
    return ipw.FloatText(**kwargs)

def FileUpload(**kwargs):
    return ipw.FileUpload(**kwargs)

def Button(**kwargs):
    return ipw.Button(**kwargs)

def HBox(items, **kwargs):
    return ipw.HBox(items)

def VBox(items, **kwargs):
    return ipw.VBox(items)

def Image(**kwargs):
    return ipw.Image(**kwargs)

def HTML(**kwargs):
    return ipw.HTML(**kwargs)

def Label(**kwargs):
    return ipw.Label(**kwargs)

def Accordion(**kwargs):
    return ipw.Accordion(**kwargs)

def Output(**kwargs):
    return ipw.Output(**kwargs)

def GridBox(items, **kwargs):
    return ipw.GridBox(items, **kwargs)

def Tab(**kwargs):
    return ipw.Tab(**kwargs)

def read_file(filename):
    return open(filename, "rb").read()

def save_file(file_content, filename):
    with open(filename, 'wb') as f:
        f.write(file_content)