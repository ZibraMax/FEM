"""Define the structure of a finite element problem. This is the parent class of individual equation classes.

The individual children classes must implement the method for calculating the element matrices and post processing.
"""


from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from .Mesh import Geometry
import os
from .Elements import *


class Core():
    def __init__(self, geometry: Geometry) -> None:
        """Create the Finite Element problem.

            Args:
                geometry (Geometry): Input geometry. The geometry must contain the elements, and the border conditions.
                You can create the geometry of the problem using the Geometry class.
        """
        self.geometry = geometry
        self.ngdl = self.geometry.ngdl
        self.K = np.zeros([self.ngdl, self.ngdl])
        self.F = np.zeros([self.ngdl, 1])
        self.Q = np.zeros([self.ngdl, 1])
        self.U = np.zeros([self.ngdl, 1])
        self.S = np.zeros([self.ngdl, 1])
        self.cbe = self.geometry.cbe
        self.cbn = self.geometry.cbn
        self.elements = self.geometry.elements

    def ensembling(self) -> None:
        """Ensembling of equation system. This method use the element gdl
        and the element matrices. The element matrices degrees of fredom must
        match the dimension of the element gdl. For m>1 variables per node,
        the gdl will be flattened. This ensure that the element matrices will always 
        be a 2-D Numpy Array.
        """
        print('Ensembling equation system...')
        for e in tqdm(self.elements, unit='Element'):
            self.K[np.ix_(e.gdlm, e.gdlm)] += e.Ke
            self.F[np.ix_(e.gdlm)] += e.Fe
            self.Q[np.ix_(e.gdlm)] += e.Qe
        print('Done!')

    def borderConditions(self) -> None:
        """Assign border conditions to the system. 
        The border conditios are assigned in this order:

        1. Natural border conditions
        2. Essential border conditions

        This ensures that in a node with 2 border conditions
        the essential border conditions will be applied.
        """
        print('Border conditions...')
        for i in tqdm(self.cbn, unit=' Natural'):
            self.Q[int(i[0])] = i[1]
        for i in tqdm(self.cbe, unit=' Essential'):
            ui = np.zeros([self.ngdl, 1])
            ui[int(i[0])] = i[1]
            vv = np.dot(self.K, ui)
            self.S = self.S - vv
            self.K[int(i[0]), :] = 0
            self.K[:, int(i[0])] = 0
            self.K[int(i[0]), int(i[0])] = 1
        self.S = self.S + self.F + self.Q
        for i in self.cbe:
            self.S[int(i[0])] = i[1]
        print('Done!')

    def solveES(self, path: str = '') -> None:
        """Solve the equation system using numpy.solve algorithm

        Args:
                path (str, optional): Path to save a text file with the solution of the problem
                This file can be loaded witouth spendign time in other finite element steps. Defaults to ''.
        """
        print('Solving equation system...')
        self.U = np.linalg.solve(self.K, self.S)
        if not path == '':
            np.savetxt(path, self.U, delimiter=',')
        for e in self.elements:
            e.setUe(self.U)
        print('Done!')

    def solve(self, path: str = '', plot: bool = True, **kargs) -> None:
        """A series of Finite Element steps

        Args:
                path (str, optional): Path to save a text file with the solution of the problem
                This file can be loaded witouth spendign time in other finite element steps. Defaults to ''.
        """
        print('Creating element matrices...')
        self.elementMatrices()
        print('Done!')
        self.ensembling()
        self.borderConditions()
        self.solveES(path)
        if plot:
            print('Post processing solution...')
            self.postProcess(**kargs)
            print('Done!')

    def solveFromFile(self, file: str, plot: bool = True, **kargs) -> None:
        """Load a solution file and show the post process for a given geometry

        Args:
                file (str): Path to the previously generated solution file.
        """
        print('Loading File...')
        self.U = np.loadtxt(file)
        for e in self.elements:
            e.setUe(self.U)
        print('Done!')
        if plot:
            print('Post processing solution...')
            self.postProcess(**kargs)
            print('Done!')

    def solveFromArray(self, solution: np.ndarray, plot: bool = True, **kargs) -> None:
        """Load a solution array to the problem.

        Args:
                solution (np.ndarray): Solution vertical array with shape (self.ngdl,1)
        """
        print('Casting solution')
        self.U = solution
        for e in self.elements:
            e.setUe(self.U)
        print('Done!')
        if plot:
            print('Post processing solution...')
            self.postProcess(**kargs)
            print('Done!')

    def profile(self) -> None:
        """Create a profile for a 3D or 2D problem.
        """
        pass

    def elementMatrices(self) -> None:
        """Calculate the element matrices
        """
        pass

    def postProcess(self) -> None:
        """Post process the solution
        """
        pass

    def print(self, filename='report.tex', aux='temp', **kargs):
        def bmatrix(a, header=[], caption='', table=False):
            rv = []
            if len(a.shape) > 2:
                raise ValueError('bmatrix can at most display two dimensions')
            if table:
                rv += [r"""\begin{table}
        \centering
        \caption{"""+caption+"""}"""]
                rv += [r'\begin{tabular}{' + '|c'*len(a[0]) + '|}\hline']
                rv += ['  ' +
                       ' & '.join([r'\textbf{'+i+'}' for i in header])+r'\\\hline']
                lines = str(a).replace('[', '').replace(
                    ']', '').splitlines()
                rv += ['  ' + ' & '.join(l.split()) +
                       r'\\ \hline' for l in lines]
                rv += [r'\end{tabular}']
                rv += [r"""\end{table}"""]
            else:
                lines = str(a).replace('[', '').replace(
                    ']', '').splitlines()
                rv += [r'\begin{bmatrix}']
                rv += ['  ' + ' & '.join(l.split()) + r'\\' for l in lines]
                rv += [r'\end{bmatrix}']
            return '\n'.join(rv)
        try:
            os.makedirs(aux)
        except Exception as e:
            print(f'{aux} folder already exists', e)
        self.geometry.show(**kargs)
        plt.savefig(aux+f'/{filename}_geometry.pdf')
        plt.close('all')
        self.solve()
        plt.savefig(aux+f'/{filename}_solution.pdf')
        plt.close('all')
        with open(f'{aux}/{filename}', 'w', encoding="utf-8") as f:
            f.write(r"""\documentclass{article}
\usepackage[utf8]{inputenc}
\usepackage{graphicx}
\usepackage{amsmath}
\usepackage{ragged2e}
\usepackage{float}
\usepackage{vmargin}
\usepackage{wrapfig}
\usepackage{titlesec}
\usepackage{lastpage}
\usepackage[spanish]{babel}
\decimalpoint

\usepackage{pgfplots}
\DeclareUnicodeCharacter{2212}{−}
\usepgfplotslibrary{groupplots,dateplot}
\usetikzlibrary{patterns,shapes.arrows}
\pgfplotsset{compat=newest}
\usepackage{hyperref}

\usepackage{caption}
\usepackage{subcaption}

\titleformat{\section}
{\sffamily}{\bfseries}{0pt}{\LARGE}

\titleformat{\subsection}
{\sffamily}{\bfseries}{0pt}{\Large}

\DeclareCaptionFont{captionFont}{\sffamily}

\captionsetup{
    width=0.9\linewidth, % width of caption is 90% of current textwidth
    labelfont=bf,        % the label, e.g. figure 12, is bold
    font=small,          % the whole caption text (label + content) is small
    labelsep=period,
    font=captionFont,
}

\setpapersize{A4}
\setmargins{3cm}        % Margen izquierda
{1cm}                   % Margen superior
{15cm}                  % Margen derecha
{22.42cm}               % Margen inferior
{10pt}
{1cm}
{0pt}
{1cm}
\usepackage{xcolor}
\usepackage{parskip}
\usepackage{fancyhdr}


\newcommand\course{\sffamily AFEM - A Python FEM Implementation}
\newcommand\hwnumber{\sffamily 2}
\newcommand\NetIDb{\sffamily \today}
\newcommand\NetIDc{\sffamily ZibraMax}
\floatplacement{figure}{H}
\floatplacement{table}{H}

\lhead{\NetIDc}
\rhead{\sffamily da.rodriguezh@uniandes.edu.co\\\course}
\lhead{\NetIDc\\\NetIDb}
	

\lfoot{\sffamily AFEM}
\cfoot{}
\rfoot{\sffamily Página \textbf{\thepage\ }de \textbf{\pageref{LastPage}}}
\headsep 1.5em

\renewcommand{\footrulewidth}{0.4pt}

\addto\captionsspanish{\renewcommand{\tablename}{Tabla}}

\begin{document}
\sffamily
\setlength{\headheight}{80pt}

\title{Reporte de resultados}
\author{}
\maketitle
\thispagestyle{fancy}
\tableofcontents
\newpage
\pagestyle{fancy}"""+'\n')
            f.write(r'\section{\textbf{Definición del modelo}}'+'\n')
            f.write(
                f'El modelo cuenta con un total del {self.ngdl} grados de libertad y {len(self.elements)} elementos. Para esta ecuación diferencial se tienen {self.geometry.nvn} variables por nodo.'+'\n')

            f.write(r'\subsection{\textbf{Elementos usados}}'+'\n')
            f.write(
                'Para este modelo se usaron los siguientes tipos de elementos: '+'\n')
            f.write(r'\begin{enumerate}'+'\n')
            nombres_elementos = {'T1V': 'Triangular de orden 1', 'T2V': 'Triangular de orden 2',
                                 'C1V': 'Cuadrilatero de orden 1', 'C2V': 'Cuadrilatero de orden 2 (Serendipity)'}
            tipos_elementos = {'T1V': LTriangular, 'T2V': QTriangular,
                               'C1V': Quadrilateral, 'C2V': Serendipity}
            for e in np.unique(self.geometry.types):
                f.write(r'\item [] \textbf{'+nombres_elementos[e]+'}\n')
                f.write(tipos_elementos[e].description()+'\n')

            f.write(r'\end{enumerate}'+'\n')

            f.write(r'\subsection{\textbf{Geometría}}'+'\n')
            f.write(r"""\begin{figure}
            \centering
            \includegraphics[width=0.95\textwidth]{"""+filename+"""_geometry.pdf}
            \caption{Geometría}
            \label{fig:figurauno}
        \end{figure}"""+'\n')

            f.write(
                r'\subsection{\textbf{Condiciones de borde escenciales}}'+'\n')
            f.write(
                r'Las condiciones de borde se dan a nivel de nodo. La numeración de los nodos se encuentra en la Figura \ref{fig:figurauno}.'+'\n')
            if len(self.cbe) > 0:
                f.write(
                    bmatrix(np.array(self.cbe), header=['Nodo', 'Valor'], caption='Condiciones de borde', table=True))
            else:
                f.write('No se definieron condiciones de borde escenciales \n')

            f.write(
                r'\subsection{\textbf{Condiciones de borde naturales}}'+'\n')
            if len(self.cbn) > 0:
                f.write(
                    bmatrix(np.array(self.cbn), header=['Nodo', 'Valor'], caption='Condiciones de borde', table=True))
            else:
                f.write('No se definieron condiciones de borde naturales \n')

            f.write(r'\section{\textbf{Solución}}'+'\n')
            f.write(r"""\begin{figure}
            \centering
            \includegraphics[width=0.95\textwidth]{"""+filename+"""_solution.pdf}
            \caption{Solución interpolada}
            \label{fig:Figurados}
        \end{figure}"""+'\n')
            f.write(
                r'\section{\textbf{Matrices y vectores de elementos}}'+'\n')
            for i, e in enumerate(self.elements):
                f.write(r'\subsection{Elemento '+format(i+1)+'}\n')
                f.write('Elemento con las siguientes coordenadas:'+'\n')
                f.write(bmatrix(e.coords, table=True,
                        header=['X', 'Y'], caption='Tabla de coordenadas para el elemento') + """\n""")
                f.write('Elemento con los siguientes grados de libertad:'+'\n')
                f.write(bmatrix(
                    e.gdl.T, table=True, header=['Variable '+format(i+1) for i in range(len(e.gdl))], caption='Tabla de grados de libertad para el elemento') + """\n""")
                f.write("""\small{$$K=""" + bmatrix(e.Ke) + """$$}\n""")
                f.write("""\small{$$F=""" + bmatrix(e.Fe) +
                        r'\hspace{0.5cm}Q='+bmatrix(e.Qe) + """$$}\n""")
                f.write("""\small{$$U_e=""" + bmatrix(e.Ue)+"""$$}\n""")
            f.write(r'\end{document}')
