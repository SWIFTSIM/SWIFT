.. Produce the header title on the second page.


.. raw:: latex

   % remove all headers/footers
   \pagestyle{empty}

   % Add background image
   \tikz[remember picture,overlay]
       \node[opacity=1.0,inner sep=0pt, anchor=north] at (current page.north)
       {\includegraphics[width=\paperwidth,height=5cm]{../../source/header_page2.jpg}};

   % Add logo
   \tikz[remember picture,overlay]
       \node[below left = -1.9cm and -1.5cm] at (current page.north)
       {\includegraphics[height=8cm]{../../source/logo/BirdW.pdf}};

   % for whatever reason, an ugly spacing is added.  Remove it here.
   \vspace{-3.3\baselineskip}
