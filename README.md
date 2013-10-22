counting
======================

This is a reference implementation of the analytical counting rules for fractional Chern insulators with an arbitrary Chern numbers.
We consider only the simplest bosonic case, corresponding to the color-entangled generalizations of Halperin 221 states.

For a given system, the code enumerates all the zero modes of the pseudopotential Hamiltonian and builds the representation of the color-entangled magnetic translation operators over these zero modes. The output is the number of zero modes in each Bloch momentum sector.

```
$ ./count -h
usage: count [-h] Ne Nx Ny C

count the number of bosonic FCI states in each momentum sector

positional arguments:
  Ne          number of bosons
  Nx          number of unit cells in the x direction
  Ny          number of unit cells in the y direction
  C           Chern number

optional arguments:
  -h, --help  show this help message and exit
```


### Examples
To reproduce the example of 2 bosons on a 3-by-2 lattice with Chern number 2 in the paper, run
```
./count 2 3 2 2
```
and the output should look like
```
ky    tot=3   
1 ├ ·  ·  ·   
0 ├ 1  1  1   
  ╰─┴──┴──┴─  
    0  1  2 kx
```
To reproduce the example of 3 bosons on a 5-by-2 lattice with Chern number 2 in the paper, run
```
./count 3 5 2 2
```
and the output should look like
```
ky       tot=10     
1 ├ 1  1  1  1  1   
0 ├ 1  1  1  1  1   
  ╰─┴──┴──┴──┴──┴─  
    0  1  2  3  4 kx
```
Finally, a non-trivial example: 6 bosons on a 6-by-5 lattice with Chern number 3:
```
$ ./count 6 6 5 3
ky          tot=2310         
4 ├ 80  75  78  76  78  75   
3 ├ 80  75  78  76  78  75   
2 ├ 80  75  78  76  78  75   
1 ├ 80  75  78  76  78  75   
0 ├ 80  75  78  76  78  75   
  ╰─┴───┴───┴───┴───┴───┴──  
    0   1   2   3   4   5  kx
```
