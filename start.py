import ipywidgets as ipw


def get_start_widget(appbase, jupbase, notebase):
    return ipw.HTML(
        f"""
    <div align="center">
        <a href="{appbase}/home.ipynb" target="_blank">
        </a>
    </div>"""
    )
