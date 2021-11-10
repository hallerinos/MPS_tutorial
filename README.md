# Welcome!

For this little tutorial, you have to implement your very own version of the iDMRG algorithm.
But don't worry, there is no need to start from scratch!
I've set up a rough code structure which is collected in a set of functions, and at its heart there are the [contractions.jl](contractions.jl) and the [iDMRG.jl](iDMRG.jl) files.
The exercises are designed such that you have to implement contractions, and add little code fragments here and there.

Before you start your intense coding session, beware to choose a unique ordering of our tensor indices beforehand, otherwise you will mix up different links in your contractions!
The choice is fully up to you, but in the solutions you will find the following one:

```
               ----                               ---- 
              |    |--(-1)                 (-1)--|    |
              |    |                             |    |
              |    |                             |    |
              |    |                             |    |
              |    |                             |    |
              |    |                             |    |
              | LW |--(-2)                 (-2)--| RW |
              |    |                             |    |
              |    |                             |    |
              |    |                             |    |
              |    |                             |    |
              |    |                             |    |
              |    |--(-3)                 (-3)--|    |
               ----                               ---- 
              
                                                  (-4)      
                                                   |        
                     -----------                  ---       
              (-1)--|    AC     |--(-4)    (-1)--| W |--(-3)
                     -----------                  ---       
                      |       |                    |        
                     (-2)    (-3)                 (-2)      
```

### The exercises are meant to be done in a certain order.
1. define the matrix product operator corresponding to the Ising Hamiltonian and check it in [isingMPO.jl](isingMPO.jl)
2. complete and check the contractions in [contractions.jl](contractions.jl)
3. scroll through [iDMRG.jl](iDMRG.jl) and substitute all occurences of whichfunhere(...) with proper functions 
4. complete the [measure.jl](measure.jl) contractions

### Here some hints
- You need to use the [TensorOperations.jl](https://github.com/Jutho/TensorOperations.jl) package.
It's enough to understand the first examples in the official [documentation](https://jutho.github.io/TensorOperations.jl/stable/indexnotation/).
- Some contractions that you know from my [TeX](notes_TeX.pdf) and [handwritten notes](notes_handwritten.pdf) are contained in [basic_contractions.jl](solutions/basic_contractions.jl).
- Have a look at [`zeros`](https://docs.julialang.org/en/v1/base/arrays/#Base.zeros) and [`rand`](https://docs.julialang.org/en/v1/stdlib/Random/#Base.rand).

In case there are any issues, do not hesitate to ask, or get some inspiration from the [solutions](solutions) folder!

Have fun! :)