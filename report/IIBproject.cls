%%%
%
%

\NeedsTeXFormat{LaTeX2e}[1995/06/01]
\ProvidesClass{IIBproject}[2000/06/05 - declaration updated. Bugs to tpl@eng]
\newif\if@cuedwide
\@cuedwidefalse

\DeclareOption{a5paper}
    {\ClassError{IIBproject}{Option `a5paper' not supported}{}}
\DeclareOption{b5paper}
    {\ClassError{IIBproject}{Option `b5paper' not supported}{}}
\DeclareOption{twocolumn}%
    {\ClassError{IIBproject}{Option `twocolumn' not supported}{}}
\DeclareOption{titlepage}%
    {\ClassError{IIBproject}{Option `titlepage' not supported}{}}
\DeclareOption{wide}
    {\@cuedwidetrue\typeout{Wide option selected}}
\PassOptionsToClass{12pt,a4paper}{report}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{report}}
\ProcessOptions
\LoadClass{report}
\if@cuedwide
\oddsidemargin=0pt \evensidemargin=0pt \topmargin=0pt \textwidth=421pt \textheight=612pt
\fi
\renewcommand{\abstractname}{Technical Abstract}
\def\projectgroup#1{\gdef\@projectgroup{#1}\gdef\theprojectgroup{#1}}
\def\name#1{\gdef\@name{#1}\gdef\thename{#1}}

\setcounter{tocdepth}{6}
\newcounter{dummy}
\newcommand\addtotoc[1]{
\refstepcounter{dummy}
\addcontentsline{toc}{chapter}{#1}}


\usepackage{graphicx}
\usepackage{setspace}
\usepackage{fullpage} % use option [cm] for a bigger effect
\usepackage[compact]{titlesec}
\titleformat{\chapter}[display]   
{\normalfont\huge\bfseries}{\chaptertitlename\ \thechapter}{10pt}{\Huge}   
\titlespacing*{\chapter}{0pt}{-50pt}{40pt}

\setlength{\headsep}{2mm}

\usepackage{fancyhdr}
\setlength{\headheight}{15.2pt}
\pagestyle{fancy}
\renewcommand{\chaptermark}[1]{ \markboth{\chaptername\ \thechapter.\ #1}{} }
\renewcommand{\sectionmark}[1]{ \markright{\thesection - #1} }

\fancypagestyle{mainstyle}{
  \fancyhf{}
  \fancyhead[L]{\textit{\nouppercase{\leftmark}} }
  \fancyhead[R]{\thepage}
}

\fancypagestyle{abstractstyle}{
  \fancyhf{}
  \fancyhead[RE]{\textit{\@title} }
  \fancyhead[LE]{\thepage}
  \fancyhead[LO]{\textit{Technical Abstract} }
  \fancyhead[RO]{\thepage}
}

\fancypagestyle{plain}{ %
  \fancyhf{} % remove everything
  \fancyfoot[C]{\thepage}
  \renewcommand{\headrulewidth}{0pt} % remove lines as well
  \renewcommand{\footrulewidth}{0pt}
}

\renewcommand\maketitle{
  \thispagestyle{empty}
  \begin{titlepage}
    \newpage
    \begin{center}
      \begin{minipage}{3.3cm}
        \makebox[\textwidth][r]{\textsc{university of}}\\
        \makebox[\textwidth][r]{\textsc{\large cambridge}}
      \end{minipage}
      \begin{minipage}[c]{2.5cm}
        \centering
        \includegraphics[height=2.5cm]{crest}
      \end{minipage}
      \begin{minipage}{3.3cm}
        \makebox[\textwidth][l]{\textsc{department of}}\\
        \makebox[\textwidth][l]{\textsc{\large engineering}}
      \end{minipage}
    \end{center}

    \vskip 4.2cm
    \begin{center}
      \begin{minipage}[t][5cm][c]{8.5cm}
        \centering
        \Large \@title \par
        \small \emph{by}\par
        \large \lineskip .5em 
        \begin{tabular}[t]{c}
          \@author
        \end{tabular}
        \par
        \small
        \emph{
          Fourth-year undergraduate project in\\Group \@projectgroup, 2012/2013
        }
      \end{minipage}
    \end{center}
    
    \vfill

    \noindent I hereby declare that, except where specifically indicated, 
    the work submitted herein is my own original work.
    \vskip 1cm %or some other suitable small distance
    \begin{center}
      Signed: \underline{\hspace{6cm}} Date: \underline{\hspace{3cm}}
    \end{center}

    \newpage
  \end{titlepage}
}

%----------------------------------------------------------------------------------------
%	ABSTRACT PAGE DESIGN
%----------------------------------------------------------------------------------------
\renewenvironment{abstract}
{
  \newpage
  \pagestyle{abstractstyle}
  \thispagestyle{plain}
  \begin{center}
    \setlength{\parskip}{0pt}
    \null
    {\large{\@title} \par}
    \smallskip
    \textit{by}\par
    \smallskip
    {\large{\@author} \par}
    \small
    \emph{
      Fourth-year undergraduate project in Group \@projectgroup,
      2012/2013}
    \par
    \bigskip
    \hrule
    \vskip 1cm
    \huge{Technical Abstract}
    \par
    \vskip 1cm
    \hrule
    \bigskip
  \end{center}
}
%----------------------------------------------------------------------------------------
{
  \vfil\vfil\vfil\null
  \cleardoublepage
}

\endinput

