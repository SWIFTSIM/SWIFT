.. Produce the header title on the first page.

.. raw:: latex

   % remove all headers/footers
   \pagestyle{empty}

   % Add image
   \tikz[remember picture,overlay]
       \node[opacity=1.0,inner sep=0pt, anchor=north] at (current page.north)
       {\includegraphics[width=\paperwidth,height=5cm]{../../source/header_page1.jpg}};

   % for whatever reason, an ugly spacing is added.  Remove it here.
   \vspace{-1.8\baselineskip}


