{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Library"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 1,
   "metadata": {},
   "outputs": [],
   "source": [
    "using Plots, Measures"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "test7 (generic function with 1 method)"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "include(\"main.jl\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "PlotLine! (generic function with 1 method)"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "function PlotLine(points::Matrix{Int64}; color::Int64=1, ls=:solid)\n",
    "    plot(points[:,1], points[:,2], points[:,3], color = color, ls=ls, lw=2,\n",
    "    zlim=(-2,2), xlim=(-2,2), ylim=(-2,2),aspect_ratio=:equal,legend= false, showaxis = false, ticks = false,\n",
    "     size = (300,300))\n",
    "end\n",
    "function PlotLine!(points::Matrix{Int64}; color::Int64=1, ls=:solid)\n",
    "    plot!(points[:,1], points[:,2], points[:,3], color = color, ls=ls, lw=2)\n",
    "end"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "background (generic function with 1 method)"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "function background()\n",
    "    points = [1 1 1; 1 -1 1; -1 -1 1; -1 1 1; 1 1 1]\n",
    "    PlotLine(points)\n",
    "    points = [1 1 -1; 1 -1 -1; -1 -1 -1]\n",
    "    PlotLine!(points)\n",
    "    points = [ -1 1 -1; 1 1 -1]\n",
    "    PlotLine!(points, ls=:dash)\n",
    "    points = [ -1 1 -1; -1 -1 -1]\n",
    "    PlotLine!(points, ls=:dash)\n",
    "    points = [1 1 -1; 1 1 1]\n",
    "    PlotLine!(points)\n",
    "    points = [1 -1 -1; 1 -1 1]\n",
    "    PlotLine!(points)\n",
    "    points = [-1 1 -1; -1 1 1]\n",
    "    PlotLine!(points, ls=:dash)\n",
    "    points = [-1 -1 -1; -1 -1 1]\n",
    "    PlotLine!(points)\n",
    "    \n",
    "    points1 = [1 1 1; 1 1 -1; 1 -1 1; 1 -1 -1; -1 1 1; -1 1 -1; -1 -1 1; -1 -1 -1]\n",
    "    points2 = [2 0 0; 0 2 0; 0 0 2; -2 0 0; 0 -2 0; 0 0 -2]\n",
    "    points = vcat(points1, points2)\n",
    "    scatter!(points[:,1], points[:,2], points[:,3], mc=:white)\n",
    "    scatter!([0], [0], [0], mc=:black)\n",
    "    end"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "PlotNeighborCoords (generic function with 1 method)"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "function PlotNeighborCoords(coords::Matrix{Int32})\n",
    "    background()\n",
    "    scatter!(coords[:,1], coords[:,2], coords[:,3])\n",
    "end"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "PlotNeighborCondition (generic function with 1 method)"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "function PlotNeighborCondition(point::Point)\n",
    "    deltas = Matrix{Int32}(undef, 0, 3)\n",
    "    for neighbor in point.neighbors\n",
    "        coord = neighbor.coord\n",
    "        delta = coord - point.coord\n",
    "        deltas = vcat(deltas, delta')\n",
    "    end\n",
    "    PlotNeighborCoords(deltas)\n",
    "end"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Test"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Initialization completed!\n"
     ]
    }
   ],
   "source": [
    "Random.seed!(1234)\n",
    "mapSize = Vector{UInt32}([300,300,300])\n",
    "universe = Universe(mapSize)\n",
    "fileName = \"/mnt/c/Users/xuke/Desktop/test6.dump\"\n",
    "RefreshFile(fileName)\n",
    "test7(universe)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 10,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "20 Vac 7 20 (152 154 148)\n"
     ]
    }
   ],
   "source": [
    "point = universe.points[20]"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 11,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "id: 5926 | probabilities: [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0] | total probability: 70.0"
     ]
    }
   ],
   "source": [
    "c = point.object.vacMigrationCondition"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 12,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "iVBORw0KGgoAAAANSUhEUgAAASwAAAEsCAIAAAD2HxkiAAAABmJLR0QA/wD/AP+gvaeTAAAbGklEQVR4nO3de1xUdf7H8c/3nIG5AYqXQcBboKaCKJKlqW3Z7trutlntlt1s01zNC3lZpYuXLHFNLbMkK3dbt6zMxyO3djd3t93WLur6KxWQQBGVvCADY4gKc2XmfH9/DCIqilw/DPN+/sHjMYjwBebF+c75zvmOkFISAPBRuAcAEOwQIQAzRAjADBECMEOEAMwQIQAzRAjADBECMEOEAMwQIQAzRAjADBECMEOEAMwQIQAzRAjADBECMEOEAMwQIQAzRAjADBECMEOEAMwQIQAzRAjADBECMEOEAMwQIQAzRAjADBECMEOEAMwQIQAzRAjADBECMEOEAMwQIQAzRAjADBECMEOEAMwQIQAzRAjADBECMNNxDwAaqbKy8p5HJ+YUWYkoqXv0x+9uCAsL4x4UNIaQUnKPARojbtiIE7983Ns/mYh0+Vk9/v524e5d3IOCxsB0NCDt3r37VPe+/gKJyNs/+VT3vrt37+YdFTQOIgxIO3fudPfoV/s97h79du7cyTUeaApEGJDGjBmjL8yr/R59Yd6YMWO4xgNNgQgDUlJSUvfy4tC9X/tvhu79unt5cVJSEu+ooHFwYiZQeb3e30yb8UVWDhHdlpz0zhuv63Q41x2QECEAM0xHAZghQgBmiBCAGSIEYIYIAZghQgBmiBCAGSIEYIYIAZghQgBmiBCAGSIEYIYIAZghQgBmiBCAGSIEYIYIAZghQgBmiBCAGSIEYIYIAZghQgBmiBCAGSIEYIYIAZghQgBmiBCAGSIEYIbX8WkTdu/efejQof79+w8dOpR7LNDacCRk5vF4Hn/88W+//TY2Nnb79u1Tpkzxer3cg4JWhZdGY7Z8+fKxY8fWHAB37dq1a9euuXPn8o4KWhOmo8yOHj1aewo6YsSIJ1/a8O+B3iijiDJStElYDBRjFlFGijKKznrGkUJLQYSc/n1S5py+6D1SynNeZU+RJKpjhhKqkMUoYkx0PlGyGEW0kaKMopuJok3CjN9nAMJ0lMfnJ+XCvb5vbDJy9zsZP+360K/G+d//4YebT1Z6b7rzIatTFtvJ/7bcI60OKnbIEkddadZiUClSTzEmEW2iyFARY6ZoY/XbSD31ChNhIa3wzUHDIMJWpUnaekJ7PlPb+4Mkoq4Gmj5Q8Xyy9GzZqX79+uXn53fv3n3BggVX+u9uH5W5Zbmb/E3635a7q1s96ZBnPfUMwKCeT/R8q7WL7WkWOpyqa3WIsJVokrYc1Z7bqx04I4koykhzEtXUBMWkIyKy2+3Hjh3r3bu3yWRqyldxei8cNi8kev7mCbus0ur5DJH66sPmRYmev9nNKBTRlAFCHRBhi/No9OERLT1bO3RWElHvcDE7QZk6QDGoDINxeqnOia7/rc1JvqveHfQqddLXMdH1JxprFh1DGz82KeXRo0c9Hk98fLxOF0SPbhFhC3L76J1D2tIsrcguiSg+QqQlKZP6KW15ylfurmOiW1Nsgx6UXv7QtLtZhF7hey8sLFy0aFFiYqLBYNizZ88TTzwxevToFvj+2iJE2CIqq+jtg9rKHK3YIYloUCcxb5DycB9FDfC5nMtHp911THT9N4vs8lxVPZ+h5kHpJRPdDWmPvvvHN8PCwohISjlx4sQ1a9Z07NixNb4rbkF00G8d5W56LU97Lc932k1EdGNXsTBZubNnO3kk5U8oxkQpVPc3dMZDVoc85aJiuyx1UqlTWp1kc8piB5U6pc1JLh8VVsjCCrpoDeZc6ZReffwFEpEQYty4cd98883YsWNb/nvihwibzSkXrf7Ot26/5j8a3NJNLExWfxLbPuq7Vh1DqWOoGEBEdVUqiWy1mix1ktUhS5301yy3PvSixRO9Xu9yuVpnzOwQYTModshVOdr6fM3hJSIa210sGKKO7hZc+V0LQeR/6s+gTkQkNEl/Oaoty9bsET2/+s93Xq+35nzM1q1bFy9ezDrY1oPHhE1yrFKu/k77w0HN6SUi+nGsSE9Rb7I0OL9l2drL3/kWDFF/N6gNn7RpPpcv2IxXMsu2ZvzszrvMJsO//vWvMbf+6P777+ceZivBkbCRCivkq7naW/ma20eKoDt7iiVD1ZQujTz6uXyy3E3OILh8okqjTUe0ZdlawVlJRL3CxJxE/4LNTfa7Enst237W7sp+ZklCj67cI209iLDBcsvlyn3apkLNq5Ei6L7rlOdTlAEdMfmsh3/BJj1LO2GXRBQXLp4arEzsp4ScP/abzeaQhNu9TurcJbieXIcIG2Dfaflyjvb+EU2TFKrQhD7KomSlbwfkVw+7l/6Yf2HBJjFSzE9qDws2zQURXpOdpfLFfb6tx6Uk0qv0m77KomSluxl3onpUVNGfDmrL9/lKnUREgzuJuYOUR/oo7WTFppkgwnrsKJErcnyfHpdEFBZCk/opTw9Wo5v0BM+g8IOLMvb7Xs3VzniIiEZGiacGt5/10uaFCK/o85Ny0V7f/9kkEUWE0LSBSlqS2gmX1dbH5qR1B3yvfFe9XurP75c9g+Ksb+Mgwkv5rzZamqXtPiWJqIuBZgxUZieqTXlqcpA4XilfvnjBZmmKOrzhCzbBBhFe4F+8WrJX21/X1UZwFd9XyDUXL9g8N1S9obELNsEG9y+iKyxeTemvGPHjqU9euVyBBZumCfZ7mUejPxdcbfEKriTntHwpR/vgiOaTFKLQhD7KwmSlHxZsGi6IIiwvL1+xYkV5eTkRxcfH/3bm7HcLdZcsXj0U36Yv9msj/lcql9dasHm8r7IwWemBBZvGCpYIfT5famrqCy+8EBcXR0Rfb98x7NGnjtyxkrB41RCXL9g8NViJMeEH1yTBEuGePXtuvfVWf4FEdMvoUUnvfFSi2Z8bHjGhr2IxEAq8Ckn0zxMyPcu3yyaJqGMopSYosxJV7IPaLIIlwqNHj9YU6Nc/vufHtqK0b69P+9anCrIYKcooYkzU1VC9o5HFSDEmYTFStKlJW6cENE3SJ8e0Zdla5vnt4WYnqjMTlIjgenZnywqWCBMTEz/99NMxY8bUvCc3N+/2B6aVeoTNJW1OsjrI6pDZZVTnrrsGtXrXXf9mu91MZDGIGDNZDKKbiaKNov2dR/VJ2lyo/T5byyuXRBRjEvOSlCn9Fewv3OyC6HrC2bNn/+hHt9511y+9Xu+6detCQkJmzpxZ869cGxz5HTknv6+g+Ai6Lpx/WuxfsPl9tnbwrCSinmFibmst2ES/X1XiJOvDId2MLf612o4gilBK+da7m1Zt2R6iU95I/fVtt9127f+35TY48r+NMoq2cEmBf3fGpVna4XOSiK4LF7MSlCcGKPrW2p0REbZ/JU6Kfr+qm5GsDzfzY5qa/Tzr3Hv32nfdvWQnz5pio02iRSP1X2206jvtpF0SUUKkSONYsAnOCDHBbx5GHcWFi7hwqnODIzo/3b18outvtdRJ5W4qd8v9Z+hKD0qv8iITPcLEtZ8p2bt376Oz55W73B30oRteejEhZfifDmov7vOVOImIkjqJ37WL3RkDCCJsJZF6itT779d13Ls9Gv3gujDRrbNV/6kjIrpSpVd5kYkeZuF/DtA333xz+7Q59ieWUMfO1nPlY2bMCxn7/Lm4W4jo5ijxNK424oAI24RQpXo/z4RIqrNSu7d6d8BSpyx2kM0prQ6qfbPWfp50eaWCyGIki1EUpc+zz1xGYRFERBGRzid/L15+ZsyLXy1MVm+LRn08EGFgMOuoT4ToE0FXmu6Wuat38iy2y1MusjpkyfkdPm1OaXNVF2uSWnWBfkazUa/89+e4G3DCT7+d6KynznoxsCPVWalPks1JNpf8+VvC4XaR3lD9D1UeM/lac5xwOTxbOSiogqJNNLiTeH3hfNP6peR2ERFVeUx/XLZsbir36IIdjoTB5e5x4/7ocqWumO9UQgw+z0tzpj/y4IPcgwp2iLBN+PdJua1Y+2msMiamxc+OPDh+/JF+9y3a65ubrExM4XiRRLgYpqNtwvYSbcU+7X+lQfTECaiBCAGYIUIAZogQgBlOzEAzKyoqWrNmTUVFRUhIyNSpUwcNGsQ9orYOEUJzslqtaWlpa9assVgsdrt93rx5kydPTklJ4R5Xm4bpKDSndevWrVq1ymKxEJHZbH711VfXr1/PPai2DhFCc7LZbLGxsTU3Q0ODdXOehkCE0Jy6dOlitVprblZV1bfjACBCaF7Tp0+fP39+WVkZETmdzrlz506aNIl7UG0dTsxAc4qNjU1PT09PT3c4HKqqTp48eejQodyDausQITSz3r17v/LKK9yjCCSYjgIwQ4QAzBAhADNE2CYkdBT3XacMjOQeB3DAiZk24YF45YF47kEAExwJAZghQgBmiBCAWUBG6HA4nlm46M7xD619/XXusQA0VeBFWFBQEHvjqFWuDluHjUvLPtEz5Sav18s9KIDGC7yzo794bPKZmcuok4WIXD37lEbFPjJ12odv/4F7XACNFHhHwnJN8Rfo50kcvisnj3E8AE0UeBEq8uLXTvC4jKHN/IqfAK0p8CIc2jM2NOd/NTfNf39n+iMBv5H7aTcVVsgzHu5xAIfAe0z46ab3ht9xZ8G2T8gSSycO333j0CenPcE9qKZ6JdeXnqUtTVEXJgfen0VoosCLUKfT7fn8Xy6Xq6CgICkpiXs4AE0VqH93DQYDCoT2IVAjBIfDkZqaarFYLBZLamqqw+HgHhE0UuBNR8EvLS3t9fNPGMrIyJBSZmRk8A4JGgdHwkC1ZcuW2jc/+ugjrpFAEyHCQCXERS8nqqp4uc9AhQgD1b333lv75j333MM1EmgiPCYMVCtXrhRCbN68mYjGjx+/YsUK7hFBIyHCQGUymdauXbt27VrugUBTYToK7UFJSUngXtGGIyEEtvfee++zzz6Li4srLi6OjY1dvHixogTYoQURQgDbuXPnsWPHNm7c6L/52WefZWRkPPnkk7yjaqgA+5sBUNuWLVtmzZpVc3Ps2LF5eYF3cSkihADmdruNRiP3KJoKEUIAS0lJ+fzzz2tunjx50mw2M46ncfCYsE14drA6N1E14rfRQI899tj06dOPHTs2YsSIw4cPb968ORCfQIsjYZtg1FGkngx45lkDKYry5ptvDhgw4KuvvlJVdePGjV26dOEeVIPhby8EvNGjR48ePZp7FI2HIyEAM0QIwAwRAjBDhADMECEAM0QIwAwRAjBDhG3Cnwu0n/zT+/5hjXsgwAAR1m/jBx/0HHZz1I2j+464ZceOHS3xJY5UyM9Pyu8rWuJzQ1uHZ8zUY/2fNsx5/xPHrJcoVG+rOHvHrKe+Wv9aSkoK97ig/cCRsB5L1q13TH6WQvVEROEd7FMWT/rdU9yDgnYFEdbDo9OTWmu+ENnlhwo733CgHUKE9dB73eStunD7tM3SIYxvONAOIcJ6LJs1w7T+BXI7iYjOng5b/8KG1au4BwXtCk7M1OOxCY+EmU3zXkxzkdIhRH3vrdeGDBnCPShoVxBh/X59772/vnjPeYBmhOkoADNECMAMEQIwQ4QAzBAhvw3vblw3eZTpzR+veXzUps2buYcDrQ1nR5m9t2nTzA2bHfNWU0iow+2anLHUaDDcPW4c97ig9eBIyGzB6rWOyQsoJJSISG9wTFk0I33VvtPS6iCf5B4ctAocCZnZSa0u0E9vOOOVQ/7iJSJFkMVAFqOIMZHFKLoZKdokLEaKMQmLkWJNokPoFT8tBBBEyOYLq0zP8jndGjntZDz/CgqV5/RCie8kbE5pc1KJk0qcMuc0EV16WHwgXtl02xW37PZJcvnIjF9vIMBvicGOEvlcpm9bsSQi0y9eML72rHPaEoqIpDNl5jeX/PONV266SUdEXo1sLlnqJKuDbE5pdVKJw39T2lwUF361L5FzWg792GvWUYxZRBkpyiiijWQximgTRRmpyC6JMN1tKxBh65FEnx7X0rO0b09JIupioBkDlVkTbj84zjxx3tNn3Z5Ig/7dP6ytuWJYp1CMScSYKLkzEYkGfa0zHjLqyO6lQ2flobN0+YGUiJZna+sPaFFGEWWkXuFiwy14KQweiLA1aJK2ntCWZGqZP0gishhp2gBl7iA1IoSIaPjw4Qd2fNm8X/G2aOF4LKSiiood0uYkq0OWOMnmlMUOsjllVhkVO6QiqMxNZW65/wzFnqvnE877xheiUDfjhQel0SbREQ9KmwMibFlVGm06ov0+Wzt4VhJRzzAxN1GZ0l9pnVdBCw+h6zuI6zvQJQfS9Cxt0V7fs0OU1IGqzSWtDvLUt8XUGwc0h/fSdxrUC+eNoo3UzUQWg4gxk8UgBncWDX1E+sUXX3TY8pHJq30SMnrqow8K0bCDf+BChC3o6xL56Je+Y5WSiPp2EE8PVib0UULa0qqQxUgWo0iMrOfDNEkZN6sn7WRzyRIHlTilzUnFDllRRccr5fFKuny6u+su3XDLFSs6affPCETNTyMjI6OqqmrfptU6ne5vf/v7nDlz1qxZ04TvLJAgwhZ0XThZHTIhUqQlKQ/FK7q2lF+DKIIm9qtj9E4vWZ2yxEE2lyw+n6jVSTanjL3qC+Yu2KO9c0gjoq4GshhFt9Aq066sv73/tv9f77nn7oKCg3l5eQkJCS3w3bQ5iLAF9TCLXXfphnQWSjudWBl1FBcu4sKpoeeNTDqKNpHNSadcdMol80qOPJN4UW/Dhg3Lzc1FhNAMhnZpp/01zbqR6rqRqibplIs+Pqq9ltUjf9vx2h9QWFg4aNAgruG1soCdIUGAk0RbT2jj/uOdttN3wGHaVy6+3l69sXJhYeGXX355ww038I6w1eBI2Hgn7HJVjhZtEs8Mxt+yBvAv2Dyfqe39QRJRVwNNH6jMeGDFhnVr3n9vIxFFRkauXbtWVYNl3RIRNkZhhXxxn/ZOgebRqJOe5iQqhmC5wzSJf8Fm+T4t/8zlCzZqWloa9wB5IMKGOXJOrszR/lSgeTVSBN3ZUzw/VEWB9fJo9OERLT1bO3RWElHvcDE7QZk6AH+8iBDhtfvutFyVo31wRPNJClFoQh9lQbJyfQecd6mH20fvHNKWZmn+J6zGR4i0JGVSvwBesGl2iLB+WWVyebb20feaJApV6KF4ZfFQpU8E8qtHZRW9fVBbkeOzOoiIBnUS8wYpD/dRVPzkLoYIr2ZHiVyR4/v0uCQis44ev15JS1JizbgT1eNcFb2xX1uZ4zvtJiJK7iyeGaL8+rr2ulzaVIiwbjtK5JJM33+LJRGFh9DEfsozQ9RuRu5htXmnXPT6ft+rudoZDxHRyCjx1GDllz0x9bwaRFiHIrsc8w9vlUad9TQrUU1NUHC5QL1KnfRKrm9tXvXzvEdGiedT1NtjcPCrHyKsQ3ezmJOodDWKJ/orYSHco2nzjlXK1d9p6/M1l48E0Z09xaJk9cauyO9aIcK6rbgR587rV1ghV+y7dMEGz9RrqOCK8PDBfPXv7zvMZttPJlksFu7hBLDccrly30ULNs8OUfp3RH6NEUSPmKfM+d3PZ82R/Ts5ulGfn929MOPP3CMKSNll8v7/+pK2eDce1lRBE/ooeb/SvXurigIbLViOhLm5uZuyDlTOSCcijahi2JiX06fPn/IYdg28dttL5LJs32dFkohMOprSX5mfpMSY0F5TBUuEGz/Y5Lz5jgu3VVUZmHzTy1/3Tbk5xiSiTeR/GxkqYszU0yzwfI7aaq+XhoXQpH7K04PVaBP3sNqLYIkwPMysnPb4ar/L4z7oDjt4XF6+NYMqyGKkqFq7p/i3JIsxi66GINrg6JLt4Trp6ckE9ckEJVLPPbL2JVginPDwwyvvm1A17Fbybx/kcoR+n/fZW8lWpyy2U7lHWh1U7Kh+a3OS1UFWh8wuozo3C9Sr1ElffdiMNla/jdRXH067mwN+b2z/1UYvZGp7al1tNDtRDZK/Pq0sWCLs1avXkkfHP798hm/QTcLjCtm/d8va1aO6iSvty1DuvtBkuZv8rdYUW+KornT/GaqzUoNKNU3Wnuj6i+1hFm1qu6faNElbjmrP7dUOnJFEFGWkOYlqaoJiCpZ7CgMhZRDtw+xyub788kuTyTRq1ChFaXwHLh+ddtdK1HFRsUV2ea6qns9gUC9O9OJio4yiRZ/l7N/ycGGysjTlwnKo/2qjZdlawVlJRL3CxJxEXG3UGoLr75vBYLjjjjvq/7h6P49avTd2yhUOpP5tyC6f6PqLPWGXLh8VVsjCCqrzQEpEkfo6Jro1xUabmrQpZ5XLQcWFVf17E0XQ+auN0rO0E3ZJRHHh4qnBuNqo9QTXkbDt8E93L5/o+lstdZJ21V9L7enu5Q9Ne4SJiCs/227ZsmV5hUXx/a4vPHyod3RXy/2LV+ZoxQ5JRImRYn4SrjZqbYiwLfJo9IOrjoluTbHl7no+Q81095KjaNY/NvXtqHvwgfH+D/tgy19n/vtU+bDfDOksnsXVRkwQYUCye8nqkCUOKnVK/2a7VgeVOqn0/KtNXGlb+35/nZ7/t/W1d5gfcd/U59e8+dNY1McmuB4TthtmHfWJEH0i6Epnd8vcVOqUpU4qtkubi6wO6U+0LFRc8hoPSZ0IBfJChO1TZz111ouBHemSSlfv779r164RI0b4b2ZmZvbu3bv1hwe1YToaXLxe7/Tp0xMSEoYMGZKbm5uZmfnGG2+EhmINnhMiDEaZmZn5+fl9+/YdNmwY91gAEQJww3IsADNECMAMEQIwQ4QAzBAhADNECMAMEQIwQ4QAzBAhADNECMAMEQIwQ4QAzBAhADNECMAMEQIwQ4QAzBAhADNECMAMEQIwQ4QAzBAhADNECMAMEQIwQ4QAzBAhADNECMAMEQIwQ4QAzBBhoPJ6vQ//dmrMDSNibhjx8G+ner1e7hFBI+FVmQLVgJE/Khx1tyflFiIK3ft13I5PDuz8intQ0Bg4EgaknJycosgYf4FE5Em5pSgyJicnh3dU0DiIMCBt27bNHZdQ+z3uuIRt27ZxjQeaAhEGpJEjR+pPFNR+j/5EwciRI7nGA02BCAPSsGHDuhYd0uVn+W/q8rO6Fh3Ca18HKJyYCVSVlZX3PDoxp8hKREndoz9+d0NYWBj3oKAxECEAM0xHAZghQgBmiBCAGSIEYIYIAZghQgBmiBCAGSIEYIYIAZghQgBmiBCAGSIEYIYIAZghQgBmiBCAGSIEYIYIAZghQgBmiBCAGSIEYIYIAZghQgBmiBCAGSIEYIYIAZghQgBmiBCAGSIEYIYIAZghQgBmiBCAGSIEYIYIAZghQgBmiBCAGSIEYIYIAZghQgBmiBCAGSIEYPb/c6vCNu6VvrAAAAAASUVORK5CYII=",
      "image/svg+xml": "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" width=\"300\" height=\"300\" viewBox=\"0 0 1200 1200\">\n<defs>\n  <clipPath id=\"clip050\">\n    <rect x=\"0\" y=\"0\" width=\"1200\" height=\"1200\"/>\n  </clipPath>\n</defs>\n<path clip-path=\"url(#clip050)\" d=\"\nM0 1200 L1200 1200 L1200 0 L0 0  Z\n  \" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\"/>\n<defs>\n  <clipPath id=\"clip051\">\n    <rect x=\"240\" y=\"120\" width=\"841\" height=\"841\"/>\n  </clipPath>\n</defs>\n<defs>\n  <clipPath id=\"clip052\">\n    <rect x=\"47\" y=\"47\" width=\"1107\" height=\"1107\"/>\n  </clipPath>\n</defs>\n<path clip-path=\"url(#clip052)\" d=\"\nM47.2441 1087.45 L47.2441 243.162 L451.89 47.2441 L1152.76 112.55 L1152.76 956.838 L748.11 1152.76 L47.2441 1087.45  Z\n  \" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\"/>\n<polyline clip-path=\"url(#clip052)\" style=\"stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:8; stroke-opacity:1; fill:none\" points=\"\n  876.378,356.275 674.055,454.234 323.622,421.581 525.945,323.622 876.378,356.275 \n  \"/>\n<polyline clip-path=\"url(#clip052)\" style=\"stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:8; stroke-opacity:1; fill:none\" points=\"\n  876.378,778.419 674.055,876.378 323.622,843.725 \n  \"/>\n<polyline clip-path=\"url(#clip052)\" style=\"stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:8; stroke-opacity:1; fill:none\" stroke-dasharray=\"32, 20\" points=\"\n  525.945,745.766 876.378,778.419 \n  \"/>\n<polyline clip-path=\"url(#clip052)\" style=\"stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:8; stroke-opacity:1; fill:none\" stroke-dasharray=\"32, 20\" points=\"\n  525.945,745.766 323.622,843.725 \n  \"/>\n<polyline clip-path=\"url(#clip052)\" style=\"stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:8; stroke-opacity:1; fill:none\" points=\"\n  876.378,778.419 876.378,356.275 \n  \"/>\n<polyline clip-path=\"url(#clip052)\" style=\"stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:8; stroke-opacity:1; fill:none\" points=\"\n  674.055,876.378 674.055,454.234 \n  \"/>\n<polyline clip-path=\"url(#clip052)\" style=\"stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:8; stroke-opacity:1; fill:none\" stroke-dasharray=\"32, 20\" points=\"\n  525.945,745.766 525.945,323.622 \n  \"/>\n<polyline clip-path=\"url(#clip052)\" style=\"stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:8; stroke-opacity:1; fill:none\" points=\"\n  323.622,843.725 323.622,421.581 \n  \"/>\n<circle clip-path=\"url(#clip052)\" cx=\"876.378\" cy=\"356.275\" r=\"14\" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n<circle clip-path=\"url(#clip052)\" cx=\"876.378\" cy=\"778.419\" r=\"14\" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n<circle clip-path=\"url(#clip052)\" cx=\"674.055\" cy=\"454.234\" r=\"14\" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n<circle clip-path=\"url(#clip052)\" cx=\"674.055\" cy=\"876.378\" r=\"14\" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n<circle clip-path=\"url(#clip052)\" cx=\"525.945\" cy=\"323.622\" r=\"14\" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n<circle clip-path=\"url(#clip052)\" cx=\"525.945\" cy=\"745.766\" r=\"14\" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n<circle clip-path=\"url(#clip052)\" cx=\"323.622\" cy=\"421.581\" r=\"14\" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n<circle clip-path=\"url(#clip052)\" cx=\"323.622\" cy=\"843.725\" r=\"14\" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n<circle clip-path=\"url(#clip052)\" cx=\"950.433\" cy=\"632.653\" r=\"14\" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n<circle clip-path=\"url(#clip052)\" cx=\"802.323\" cy=\"502.041\" r=\"14\" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n<circle clip-path=\"url(#clip052)\" cx=\"600\" cy=\"177.856\" r=\"14\" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n<circle clip-path=\"url(#clip052)\" cx=\"249.567\" cy=\"567.347\" r=\"14\" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n<circle clip-path=\"url(#clip052)\" cx=\"397.677\" cy=\"697.959\" r=\"14\" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n<circle clip-path=\"url(#clip052)\" cx=\"600\" cy=\"1022.14\" r=\"14\" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n<circle clip-path=\"url(#clip052)\" cx=\"600\" cy=\"600\" r=\"14\" fill=\"#000000\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"none\"/>\n<circle clip-path=\"url(#clip052)\" cx=\"249.567\" cy=\"567.347\" r=\"14\" fill=\"#00a8cb\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n<circle clip-path=\"url(#clip052)\" cx=\"397.677\" cy=\"697.959\" r=\"14\" fill=\"#00a8cb\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n<circle clip-path=\"url(#clip052)\" cx=\"323.622\" cy=\"843.725\" r=\"14\" fill=\"#00a8cb\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n<circle clip-path=\"url(#clip052)\" cx=\"674.055\" cy=\"454.234\" r=\"14\" fill=\"#00a8cb\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n<circle clip-path=\"url(#clip052)\" cx=\"600\" cy=\"1022.14\" r=\"14\" fill=\"#00a8cb\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n<circle clip-path=\"url(#clip052)\" cx=\"525.945\" cy=\"745.766\" r=\"14\" fill=\"#00a8cb\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n<circle clip-path=\"url(#clip052)\" cx=\"600\" cy=\"177.856\" r=\"14\" fill=\"#00a8cb\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n</svg>\n",
      "text/html": [
       "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n",
       "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" width=\"300\" height=\"300\" viewBox=\"0 0 1200 1200\">\n",
       "<defs>\n",
       "  <clipPath id=\"clip080\">\n",
       "    <rect x=\"0\" y=\"0\" width=\"1200\" height=\"1200\"/>\n",
       "  </clipPath>\n",
       "</defs>\n",
       "<path clip-path=\"url(#clip080)\" d=\"\n",
       "M0 1200 L1200 1200 L1200 0 L0 0  Z\n",
       "  \" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\"/>\n",
       "<defs>\n",
       "  <clipPath id=\"clip081\">\n",
       "    <rect x=\"240\" y=\"120\" width=\"841\" height=\"841\"/>\n",
       "  </clipPath>\n",
       "</defs>\n",
       "<defs>\n",
       "  <clipPath id=\"clip082\">\n",
       "    <rect x=\"47\" y=\"47\" width=\"1107\" height=\"1107\"/>\n",
       "  </clipPath>\n",
       "</defs>\n",
       "<path clip-path=\"url(#clip082)\" d=\"\n",
       "M47.2441 1087.45 L47.2441 243.162 L451.89 47.2441 L1152.76 112.55 L1152.76 956.838 L748.11 1152.76 L47.2441 1087.45  Z\n",
       "  \" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\"/>\n",
       "<polyline clip-path=\"url(#clip082)\" style=\"stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:8; stroke-opacity:1; fill:none\" points=\"\n",
       "  876.378,356.275 674.055,454.234 323.622,421.581 525.945,323.622 876.378,356.275 \n",
       "  \"/>\n",
       "<polyline clip-path=\"url(#clip082)\" style=\"stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:8; stroke-opacity:1; fill:none\" points=\"\n",
       "  876.378,778.419 674.055,876.378 323.622,843.725 \n",
       "  \"/>\n",
       "<polyline clip-path=\"url(#clip082)\" style=\"stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:8; stroke-opacity:1; fill:none\" stroke-dasharray=\"32, 20\" points=\"\n",
       "  525.945,745.766 876.378,778.419 \n",
       "  \"/>\n",
       "<polyline clip-path=\"url(#clip082)\" style=\"stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:8; stroke-opacity:1; fill:none\" stroke-dasharray=\"32, 20\" points=\"\n",
       "  525.945,745.766 323.622,843.725 \n",
       "  \"/>\n",
       "<polyline clip-path=\"url(#clip082)\" style=\"stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:8; stroke-opacity:1; fill:none\" points=\"\n",
       "  876.378,778.419 876.378,356.275 \n",
       "  \"/>\n",
       "<polyline clip-path=\"url(#clip082)\" style=\"stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:8; stroke-opacity:1; fill:none\" points=\"\n",
       "  674.055,876.378 674.055,454.234 \n",
       "  \"/>\n",
       "<polyline clip-path=\"url(#clip082)\" style=\"stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:8; stroke-opacity:1; fill:none\" stroke-dasharray=\"32, 20\" points=\"\n",
       "  525.945,745.766 525.945,323.622 \n",
       "  \"/>\n",
       "<polyline clip-path=\"url(#clip082)\" style=\"stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:8; stroke-opacity:1; fill:none\" points=\"\n",
       "  323.622,843.725 323.622,421.581 \n",
       "  \"/>\n",
       "<circle clip-path=\"url(#clip082)\" cx=\"876.378\" cy=\"356.275\" r=\"14\" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n",
       "<circle clip-path=\"url(#clip082)\" cx=\"876.378\" cy=\"778.419\" r=\"14\" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n",
       "<circle clip-path=\"url(#clip082)\" cx=\"674.055\" cy=\"454.234\" r=\"14\" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n",
       "<circle clip-path=\"url(#clip082)\" cx=\"674.055\" cy=\"876.378\" r=\"14\" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n",
       "<circle clip-path=\"url(#clip082)\" cx=\"525.945\" cy=\"323.622\" r=\"14\" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n",
       "<circle clip-path=\"url(#clip082)\" cx=\"525.945\" cy=\"745.766\" r=\"14\" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n",
       "<circle clip-path=\"url(#clip082)\" cx=\"323.622\" cy=\"421.581\" r=\"14\" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n",
       "<circle clip-path=\"url(#clip082)\" cx=\"323.622\" cy=\"843.725\" r=\"14\" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n",
       "<circle clip-path=\"url(#clip082)\" cx=\"950.433\" cy=\"632.653\" r=\"14\" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n",
       "<circle clip-path=\"url(#clip082)\" cx=\"802.323\" cy=\"502.041\" r=\"14\" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n",
       "<circle clip-path=\"url(#clip082)\" cx=\"600\" cy=\"177.856\" r=\"14\" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n",
       "<circle clip-path=\"url(#clip082)\" cx=\"249.567\" cy=\"567.347\" r=\"14\" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n",
       "<circle clip-path=\"url(#clip082)\" cx=\"397.677\" cy=\"697.959\" r=\"14\" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n",
       "<circle clip-path=\"url(#clip082)\" cx=\"600\" cy=\"1022.14\" r=\"14\" fill=\"#ffffff\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n",
       "<circle clip-path=\"url(#clip082)\" cx=\"600\" cy=\"600\" r=\"14\" fill=\"#000000\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"none\"/>\n",
       "<circle clip-path=\"url(#clip082)\" cx=\"249.567\" cy=\"567.347\" r=\"14\" fill=\"#00a8cb\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n",
       "<circle clip-path=\"url(#clip082)\" cx=\"397.677\" cy=\"697.959\" r=\"14\" fill=\"#00a8cb\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n",
       "<circle clip-path=\"url(#clip082)\" cx=\"323.622\" cy=\"843.725\" r=\"14\" fill=\"#00a8cb\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n",
       "<circle clip-path=\"url(#clip082)\" cx=\"674.055\" cy=\"454.234\" r=\"14\" fill=\"#00a8cb\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n",
       "<circle clip-path=\"url(#clip082)\" cx=\"600\" cy=\"1022.14\" r=\"14\" fill=\"#00a8cb\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n",
       "<circle clip-path=\"url(#clip082)\" cx=\"525.945\" cy=\"745.766\" r=\"14\" fill=\"#00a8cb\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n",
       "<circle clip-path=\"url(#clip082)\" cx=\"600\" cy=\"177.856\" r=\"14\" fill=\"#00a8cb\" fill-rule=\"evenodd\" fill-opacity=\"1\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-width=\"2.4\"/>\n",
       "</svg>\n"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "PlotNeighborCondition(point)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Julia 1.7.2",
   "language": "julia",
   "name": "julia-1.7"
  },
  "language_info": {
   "file_extension": ".jl",
   "mimetype": "application/julia",
   "name": "julia",
   "version": "1.7.2"
  },
  "orig_nbformat": 4
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
