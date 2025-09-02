import os
import sys
from datetime import datetime

project = "STRmie-HD"
author = "Mazzalab"
year = str(datetime.now().year)
release = "0.1.0"

extensions = [
    "myst_parser",
    "sphinx_copybutton",
]

myst_enable_extensions = [
    "colon_fence",
    "linkify",
    "deflist",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "sphinx_rtd_theme"
html_theme_options = {
    "navigation_depth": 3,
    "collapse_navigation": False,
}
html_static_path = []
html_logo = None
html_title = f"{project} Documentation"
