"""Setting up CP2K plugin for AiiDA"""

import json

from setuptools import find_packages, setup


def run_setup():
    with open("setup.json", encoding="utf-8") as info:
        kwargs = json.load(info)
    setup(
        include_package_data=True,
        packages=find_packages(),
        long_description=open("README.md", encoding="utf-8").read(),
        long_description_content_type="text/markdown",
        **kwargs,
    )


if __name__ == "__main__":
    run_setup()
