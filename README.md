# Welcome!

Before you start digging in the exercises, let's fix some necessities right away.
First of all, we must choose a unique ordering of our tensors that we want to contract in a network.
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

The exercises are meant to be done in a certain order
1. define the matrix product operator in [isingMPO.jl](isingMPO.jl)
2. complete the contractions in [contractions.jl](contractions.jl)
3. scroll through [iDMRG.jl](iDMRG.jl) and substitute all occurences of whichfunhere(...) with proper functions 
4. complete the [measure.jl](measure.jl) contractions

In case there are any issues, do not hesitate to ask, or get some inspiration from the [solutions](solutions) folder!