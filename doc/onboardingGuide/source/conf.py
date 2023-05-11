# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = "SWIFT Onboarding Guide"
copyright = "2023, SWIFT Collaboration"
author = "SWIFT Collaboration"

# The full version, including alpha/beta/rc tags
release = "1.0"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = []

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "alabaster"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]


# -- Options for LaTeX-PDF output --------------------------------------------
# See https://www.sphinx-doc.org/en/master/latex.html for more options

latex_theme = "howto"

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    "papersize": "a4paper",
    # The font size ('10pt', '11pt' or '12pt').
    "pointsize": "10pt",
    # Font package inclusion.
    "fontpkg": r"""
\usepackage[sfdefault]{AlegreyaSans} %% Option 'black' gives heavier bold face
\renewcommand*\oldstylenums[1]{{\AlegreyaSansOsF #1}}
""",
    # Other possible choices:
    # \usepackage[sfdefault]{ClearSans} %% option 'sfdefault' activates Clear Sans as the default text font
    # \usepackage[sfdefault,light]{merriweather} %% Option 'black' gives heavier bold face
    #  %\usepackage[sfdefault,book]{FiraSans} %% option 'sfdefault' activates Fira Sans as the default text font
    #  %\usepackage{marcellus} %% option 'sfdefault' activates Fira Sans as the default text font
    # ------------------------------------------
    # override use of fncychap
    "fncychap": r"""
""",
    # pass options to packages sphinx already loads
    "passoptionstopackages": r"""
\PassOptionsToPackage{top=6.1cm, bottom=1cm, left=0.6cm, right=0.6cm}{geometry}
""",
    # Additional stuff for the LaTeX preamble.
    "preamble": r"""
\usepackage{multicol} % use two columns throughout the document

\setcounter{secnumdepth}{0} % turn off chapter numbering

% Reduce spacing after headings
%------------------------------------
% https://tex.stackexchange.com/questions/53338/reducing-spacing-after-headings
\usepackage{titlesec}

\titlespacing\title{0pt}{0pt plus 0pt minus 0pt}{0pt plus 0pt minus 0pt}
\titlespacing\chapter{0pt}{0pt plus 0pt minus 0pt}{0pt plus 0pt minus 0pt}
\titlespacing\section{0pt}{4pt plus 0pt minus 8pt}{4pt plus 0pt minus 4pt}
\titlespacing\subsection{0pt}{1pt plus 0pt minus 4pt}{0pt plus 0pt minus 8pt}
\titlespacing\subsubsection{0pt}{1pt plus 0pt minus 2pt}{0pt plus 0pt minus 8pt}

% Reduce section font sizes
%------------------------------------
\titleformat*{\section}{\Large\bf}
\titleformat*{\subsection}{\normalsize\bf}
\titleformat*{\subsubsection}{\small\bf}



% Modify the way inline verbatim behaves.
%------------------------------------------
% this version changes the text color.
%\definecolor{inlineVerbatimTextColor}{rgb}{0.6, 0.4, 0.5}
%\protect\renewcommand{\sphinxcode}[1]{\textcolor{inlineVerbatimTextColor}{\texttt{#1}}}

% this version keeps the text color, but adds a colorful box.
\definecolor{inlineVerbatimBorderColor}{rgb}{0.90, 1.00, 0.95}

\protect\renewcommand{\sphinxcode}[1]{\colorbox{inlineVerbatimBorderColor}{\texttt{#1}}}


% Drawing and image positioning
%---------------------------------
\usepackage{tikz}
\usetikzlibrary{positioning}


% Reduce space between \item
%----------------------------------

\let\tempone\itemize
\let\temptwo\enditemize
\renewenvironment{itemize}{\tempone\addtolength{\itemsep}{-0.5\baselineskip}}{\temptwo}
""",
    #  last thing before \begin{document}:
    "makeindex": r"""
\pagestyle{empty}
""",
    #  override title making.
    "maketitle": r"""
\begin{multicols}{2} % make two columns
""",
    #  additional footer content
    "atendofbody": r"""
\end{multicols} % make two columns
""",
    #  override ToC.
    "tableofcontents": r"""
""",
    "extraclassoptions": r"""
""",
    #  sphinx related stuff
    "sphinxsetup": r"""
VerbatimColor={rgb}{0.90, 1.00, 0.95},
verbatimwithframe=false,
hmargin={0.6cm, 0.6cm},
vmargin={5.2cm, 1cm},
TitleColor={rgb}{0.40,0.00,0.33},
OuterLinkColor={rgb}{0., 0.40, 0.27},
""",
}

latex_documents = [
    (
        "index",  # startdocname
        "onboarding.tex",  # targetname
        "SWIFT Onboarding Guide",  # title
        "SWIFT Collaboration",  # author
        "howto",  # theme
        False,  # toctree only
    )
]
